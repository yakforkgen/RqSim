# GPU Performance Analysis Report

## Проблема
При 500 узлах GPU загружен минимально (~5-10%), скорость симуляции низкая.

## Выявленные узкие места

### ?? КРИТИЧЕСКИЕ ПРОБЛЕМЫ

#### 1. Избыточные CPU-GPU трансферы данных (Form_Main.cs, строки 740-780)

**Каждый шаг** выполняется:
```csharp
// Строки 745-760: КАЖДЫЕ 10 шагов
float[] hostWeights = new float[edgeCount];
float[] hostMasses = new float[graph.N];
for (int e = 0; e < edgeCount; e++) { ... }  // O(E) на CPU!
_gpuGravityEngine.UploadInitialData(...);    // CPU ? GPU transfer
```

**Проблема**: При 500 узлах и ~18000 рёбер:
- 18000 ? sizeof(float) = 72 KB на каждый upload
- Аллокация new float[] создаёт GC pressure
- CPU loop O(E) выполняется синхронно

**Решение**: Выделить буферы один раз, переиспользовать

#### 2. Синхронизация weights обратно на CPU (строки 768-780)

```csharp
// КАЖДЫЙ VisualizationInterval (1 шаг по умолчанию!)
float[] hostWeights = new float[edgeCount];
_gpuGravityEngine.SyncToHost(hostWeights);  // GPU ? CPU
for (int e = 0; e < edgeCount; e++) { ... } // O(E) на CPU
```

**Проблема**: 
- SyncToHost вызывает GPU sync barrier ? GPU простаивает
- VisualizationInterval=1 означает sync КАЖДЫЙ шаг

**Решение**: Увеличить VisualizationInterval до 10-50

#### 3. Топология перестраивается каждые 10 шагов (строки 792-815)

```csharp
if (step % 10 == 0) {
    graph.Step();
    graph.BuildSoAViews();              // Перестройка CSR на CPU
    _gpuGravityEngine.UpdateTopologyBuffers(graph);  // Re-upload
    
    // Ещё хуже - для ScalarFieldEngine:
    for (int idx = 0; idx < totalDirectedEdges; idx++) {
        // ВЛОЖЕННЫЙ ЦИКЛ O(N) для поиска from!
        for (int n = 0; n < graph.N; n++) {
            if (idx >= graph.CsrOffsets[n] && idx < graph.CsrOffsets[n + 1]) {
                from = n; break;
            }
        }
    }
}
```

**Проблема**: O(E ? N) = O(18000 ? 500) = 9M операций!

**Решение**: Предвычислить mapping CSR index ? node

#### 4. ScalarFieldEngine: копирование каждый шаг (строки 830-845)

```csharp
float[] scalarField = new float[graph.N];   // Alloc!
for (int n = 0; n < graph.N; n++)
    scalarField[n] = (float)graph.ScalarField[n];  // Copy!

_gpuScalarFieldEngine.UpdateField(...);

for (int n = 0; n < graph.N; n++)
    graph.ScalarField[n] = scalarField[n];  // Copy back!
```

**Проблема**: 2 ? N копирований + alloc каждый шаг

**Решение**: Использовать UpdateFieldNoCopy + SyncToHost периодически

#### 5. Статистика каждый шаг (строки 860-875)

```csharp
int excited = graph.State.Count(s => s == NodeState.Excited); // LINQ!
var heavyStats = graph.GetHeavyClusterStatsCorrelationMass(...);
var clusters = graph.GetStrongCorrelationClusters(...);
```

**Проблема**: Тяжёлые CPU операции каждый шаг, GPU ждёт

### ?? СРЕДНИЕ ПРОБЛЕМЫ

#### 6. Маленький размер работы для GPU

При 500 узлах:
- ~18000 рёбер ? 18000 threads для GravityShader
- ThreadGroupSize = 64 ? 282 thread groups
- GTX 1650 имеет 14 SM ? 64 CUDA cores = 896 cores

**Проблема**: GPU недозагружен, но overhead запуска ядра фиксирован

#### 7. SpectralWalkEngine: 50000 walkers но только 100 шагов

```csharp
int[] returns = _gpuSpectralWalkEngine.RunSteps(100);
```

**Проблема**: Каждый Step() делает sync ? 100 syncs для spectral dimension

### ?? РЕКОМЕНДАЦИИ

#### Быстрые исправления (2-4x ускорение)

1. **Увеличить VisualizationInterval до 10-50**
2. **Выделить буферы один раз** в начале симуляции
3. **Предвычислить CSR node mapping** вместо вложенного цикла
4. **Уменьшить частоту topology rebuild** (каждые 50-100 шагов)

#### Архитектурные улучшения (10x+ ускорение)

1. **Fused kernel**: Объединить curvature + gravity + weight update в один шейдер
2. **Persistent GPU state**: Держать все данные на GPU, sync только для визуализации
3. **Batched statistics**: Вычислять метрики на GPU (excited count, cluster stats)
4. **Async copy**: Использовать async transfers где возможно

## Пример оптимизированного кода

```csharp
// Предвыделенные буферы
private float[]? _hostWeightsBuffer;
private float[]? _hostMassesBuffer;
private int[]? _csrNodeMapping;  // CSR index ? source node

// Инициализация (один раз)
_hostWeightsBuffer = new float[edgeCount];
_hostMassesBuffer = new float[graph.N];
_csrNodeMapping = BuildCsrNodeMapping(graph);

// Главный цикл
for (int step = 0; step < totalSteps; step++)
{
    // GPU step без sync
    _gpuGravityEngine.EvolveFullGpuStep_NoCopy(...);
    
    // Sync только для визуализации
    if (step % 50 == 0)
    {
        _gpuGravityEngine.SyncToHost(_hostWeightsBuffer);
        // Update graph.Weights
    }
    
    // Topology rebuild редко
    if (step % 100 == 0 && topologyChanged)
    {
        UpdateTopologyOptimized();
    }
}
```

## Ожидаемое улучшение

| Оптимизация | Ускорение |
|-------------|-----------|
| Buffer reuse | 1.5x |
| Reduce sync frequency | 2-3x |
| CSR node mapping | 2x |
| Fused kernels | 2-3x |
| **Итого** | **10-20x** |

Для 500 узлов при текущей реализации ~1 step/sec.
После оптимизации ожидается 10-20 steps/sec.
