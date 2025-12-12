using System;
using System.Collections.Generic;

namespace RQSimulation
{
    /// <summary>
    /// Event types for the event-based simulation
    /// </summary>
    public enum NodeEventType
    {
        Update,      // Node state update
        Signal,      // Signal arrival from neighbor
        Measurement  // Measurement event
    }

    /// <summary>
    /// Event structure for priority queue-based simulation.
    /// Implements RQ-hypothesis checklist item 2: True Relational Time.
    /// </summary>
    public readonly struct NodeUpdateEvent
    {
        public readonly double Time;      // Global reference time when event occurs
        public readonly int NodeId;       // Node to update
        public readonly NodeEventType ActionType;
        public readonly int SourceNode;   // Source node for signals (-1 if not applicable)

        public NodeUpdateEvent(double time, int nodeId, NodeEventType actionType, int sourceNode = -1)
        {
            Time = time;
            NodeId = nodeId;
            ActionType = actionType;
            SourceNode = sourceNode;
        }
    }

    public partial class RQGraph
    {
        // Node-specific time tracking for asynchronous updates
        private double[]? _nodeProperTime;
        private double[]? _nodeNextUpdateTime;
        private double[]? _nodeTimeDilationFactor;
        private double _globalTime = 0.0;

        // Event-based simulation queue (checklist item 2)
        private PriorityQueue<NodeUpdateEvent, double>? _eventQueue;

        /// <summary>
        /// Current reference time for event-based simulation.
        /// Implements checklist item 2.2: Global coordinate time.
        /// </summary>
        public double CurrentReferenceTime
        {
            get => _globalTime;
            set => _globalTime = value;
        }

        // Configuration
        private const double BaseTimestep = 0.01;
        private const double MinTimeDilation = 0.1;  // Minimum time flow rate
        private const double MaxTimeDilation = 1.0;  // Maximum (normal) time flow

        /// <summary>
        /// Initialize asynchronous time tracking for nodes
        /// </summary>
        public void InitAsynchronousTime()
        {
            _nodeProperTime = new double[N];
            _nodeNextUpdateTime = new double[N];
            _nodeTimeDilationFactor = new double[N];
            _eventQueue = new PriorityQueue<NodeUpdateEvent, double>();

            // Initialize all nodes to start together
            for (int i = 0; i < N; i++)
            {
                _nodeProperTime[i] = 0.0;
                _nodeNextUpdateTime[i] = 0.0;
                _nodeTimeDilationFactor[i] = 1.0; // Normal time flow initially

                // Schedule initial update events for all nodes
                _eventQueue.Enqueue(new NodeUpdateEvent(0.0, i, NodeEventType.Update), 0.0);
            }

            _globalTime = 0.0;
        }

        /// <summary>
        /// Schedule an event for a node at a specific time.
        /// Implements checklist item 2.3: Event scheduling.
        /// </summary>
        /// <param name="nodeId">Node to schedule event for</param>
        /// <param name="time">Time at which event should occur</param>
        /// <param name="eventType">Type of event</param>
        /// <param name="sourceNode">Source node for signal events</param>
        public void Schedule(int nodeId, double time, NodeEventType eventType = NodeEventType.Update, int sourceNode = -1)
        {
            if (_eventQueue == null)
                InitAsynchronousTime();

            var evt = new NodeUpdateEvent(time, nodeId, eventType, sourceNode);
            _eventQueue!.Enqueue(evt, time);
        }

        /// <summary>
        /// Propagate a signal from source to target with finite speed of light delay.
        /// Implements checklist item 2.4: Causal signal propagation.
        /// </summary>
        /// <param name="source">Source node</param>
        /// <param name="target">Target node</param>
        public void PropagateSignal(int source, int target)
        {
            if (_eventQueue == null)
                InitAsynchronousTime();

            // Validate speed of light is positive to prevent division by zero
            double c = PhysicsConstants.SpeedOfLight;
            if (c <= 0)
                c = 1.0; // Fallback to default

            // Calculate arrival time based on graph distance
            double dist = GetGraphDistanceWeighted(source, target);
            double arrivalTime = CurrentReferenceTime + (dist / c);

            // Schedule signal arrival event
            Schedule(target, arrivalTime, NodeEventType.Signal, source);
        }

        /// <summary>
        /// Process a single node event.
        /// Implements checklist item 2.5: Event processing.
        /// </summary>
        /// <param name="nodeId">Node to process</param>
        /// <param name="eventType">Type of event</param>
        /// <param name="sourceNode">Source node for signal events</param>
        public void ProcessNodeEvent(int nodeId, NodeEventType eventType = NodeEventType.Update, int sourceNode = -1)
        {
            if (nodeId < 0 || nodeId >= N)
                return;

            switch (eventType)
            {
                case NodeEventType.Update:
                    // Regular node update with proper time
                    double dTau = BaseTimestep * GetNodeTimeDilation(nodeId);
                    UpdateSingleNode(nodeId, dTau);

                    // Advance proper time
                    if (_nodeProperTime != null)
                        _nodeProperTime[nodeId] += dTau;

                    // Schedule next update
                    Schedule(nodeId, CurrentReferenceTime + BaseTimestep, NodeEventType.Update);
                    break;

                case NodeEventType.Signal:
                    // Handle signal arrival from another node
                    ProcessSignalArrival(nodeId, sourceNode);
                    break;

                case NodeEventType.Measurement:
                    // Handle measurement event
                    ProcessMeasurementEvent(nodeId);
                    break;
            }
        }

        /// <summary>
        /// Process signal arrival at a node.
        /// </summary>
        private void ProcessSignalArrival(int target, int source)
        {
            if (source < 0 || source >= N || target < 0 || target >= N)
                return;

            // Signal can trigger excitation if source was excited
            if (State[source] == NodeState.Excited || State[source] == NodeState.Refractory)
            {
                // Probability of excitation based on connection strength
                double connectionStrength = Edges[source, target] ? Weights[source, target] : 0;
                if (connectionStrength > 0 && _rng.NextDouble() < connectionStrength * PhysicsConstants.SignalExcitationProbability)
                {
                    if (State[target] == NodeState.Rest)
                    {
                        State[target] = NodeState.Excited;

                        // Propagate to neighbors (causal cone expands)
                        foreach (int neighbor in Neighbors(target))
                        {
                            PropagateSignal(target, neighbor);
                        }
                    }
                }
            }

            // Update local potential based on signal
            if (LocalPotential != null && target < LocalPotential.Length)
            {
                double signalStrength = Edges[source, target] ? Weights[source, target] * PhysicsConstants.SignalStrengthFactor : 0;
                LocalPotential[target] += signalStrength;
            }
        }

        /// <summary>
        /// Process measurement event at a node.
        /// </summary>
        private void ProcessMeasurementEvent(int nodeId)
        {
            // Collapse quantum state at node
            if (_waveMulti != null)
            {
                int d = GaugeDimension;
                int baseIdx = nodeId * d;

                // Measure and collapse
                double totalProb = 0;
                for (int a = 0; a < d && baseIdx + a < _waveMulti.Length; a++)
                {
                    totalProb += _waveMulti[baseIdx + a].Magnitude * _waveMulti[baseIdx + a].Magnitude;
                }

                // Normalize to collapse
                if (totalProb > 1e-10)
                {
                    double invNorm = 1.0 / Math.Sqrt(totalProb);
                    for (int a = 0; a < d && baseIdx + a < _waveMulti.Length; a++)
                    {
                        _waveMulti[baseIdx + a] *= invNorm;
                    }
                }
            }
        }

        /// <summary>
        /// Run event-based simulation for specified duration.
        /// Implements checklist item 2.1: Event queue processing.
        /// </summary>
        /// <param name="duration">How long to simulate</param>
        /// <param name="maxEvents">Maximum events to process (for safety)</param>
        public void RunEventBased(double duration, int maxEvents = 100000)
        {
            if (_eventQueue == null)
                InitAsynchronousTime();

            double endTime = CurrentReferenceTime + duration;
            int eventCount = 0;

            while (_eventQueue!.Count > 0 && eventCount < maxEvents)
            {
                // Peek at next event
                if (!_eventQueue.TryPeek(out var evt, out double time))
                    break;

                // Check if we've exceeded duration
                if (time > endTime)
                    break;

                // Dequeue and process
                _eventQueue.Dequeue();
                CurrentReferenceTime = time;
                ProcessNodeEvent(evt.NodeId, evt.ActionType, evt.SourceNode);

                eventCount++;
            }
        }

        /// <summary>
        /// Run relational dynamics loop with Page-Wootters style time evolution.
        /// Each node updates at a rate determined by its local curvature: dτ = 1/√|R| + ε
        /// Higher curvature (gravity wells) results in slower local time (time dilation).
        /// Implements checklist item 2.1: Page-Wootters dynamics with curvature-based scheduling.
        /// </summary>
        /// <param name="maxEvents">Maximum events to process</param>
        public void RunRelationalLoop(int maxEvents = 100000)
        {
            if (_eventQueue == null)
                InitAsynchronousTime();

            int eventCount = 0;

            while (_eventQueue!.TryDequeue(out var evt, out double time) && eventCount < maxEvents)
            {
                // "Now" is the time of the current event
                CurrentReferenceTime = time;
                int nodeId = evt.NodeId;

                // Process node update
                UpdateLocalNode(nodeId);

                // Schedule next update based on local curvature (gravity = time dilation)
                // dTau = 1 / sqrt(|R_local| + epsilon) - higher curvature = slower time
                double localCurvature = GetLocalCurvature(nodeId);
                double dTau = 1.0 / Math.Sqrt(Math.Abs(localCurvature) + PhysicsConstants.CurvatureRegularizationEpsilon);

                // Clamp dTau to reasonable range
                dTau = Math.Clamp(dTau, BaseTimestep * MinTimeDilation, BaseTimestep / MinTimeDilation);

                // Schedule next event for this node
                _eventQueue.Enqueue(new NodeUpdateEvent(time + dTau, nodeId, NodeEventType.Update), time + dTau);

                eventCount++;
            }
        }

        /// <summary>
        /// Update a single node's state during relational loop.
        /// Called by RunRelationalLoop for each scheduled event.
        /// </summary>
        private void UpdateLocalNode(int nodeId)
        {
            if (nodeId < 0 || nodeId >= N)
                return;

            // Advance proper time for this node
            if (_nodeProperTime != null)
            {
                double dTau = BaseTimestep * GetNodeTimeDilation(nodeId);
                _nodeProperTime[nodeId] += dTau;
            }

            // Update excitation state
            if (State[nodeId] == NodeState.Refractory)
            {
                _refractoryCounter[nodeId]--;
                if (_refractoryCounter[nodeId] <= 0)
                    State[nodeId] = NodeState.Rest;
                return;
            }

            // Check for neighbor-driven excitation
            double excitationStrength = 0.0;
            foreach (int j in Neighbors(nodeId))
            {
                if (State[j] == NodeState.Excited)
                {
                    excitationStrength += Weights[nodeId, j];
                }
            }

            // Probabilistic excitation based on neighbor activity
            if (State[nodeId] == NodeState.Rest && excitationStrength > 0.5)
            {
                if (_rng.NextDouble() < excitationStrength)
                {
                    State[nodeId] = NodeState.Excited;

                    // Propagate signals to neighbors (causal propagation)
                    foreach (int j in Neighbors(nodeId))
                    {
                        PropagateSignal(nodeId, j);
                    }
                }
            }
            else if (State[nodeId] == NodeState.Excited)
            {
                // Transition to refractory
                State[nodeId] = NodeState.Refractory;
                _refractoryCounter[nodeId] = DynamicBaseRefractorySteps;
            }

            // Update local potential field
            if (LocalPotential != null && nodeId < LocalPotential.Length)
            {
                // Diffusion from neighbors
                double diffusion = 0.0;
                int neighborCount = 0;
                foreach (int j in Neighbors(nodeId))
                {
                    if (j < LocalPotential.Length)
                    {
                        diffusion += (LocalPotential[j] - LocalPotential[nodeId]) * Weights[nodeId, j];
                        neighborCount++;
                    }
                }
                if (neighborCount > 0)
                {
                    double dTau = BaseTimestep * GetNodeTimeDilation(nodeId);
                    LocalPotential[nodeId] += 0.1 * dTau * diffusion / neighborCount;
                }

                // Natural decay
                LocalPotential[nodeId] *= 0.99;
            }
        }

        /// <summary>
        /// Update time dilation factors based on local mass and curvature
        /// </summary>
        public void UpdateTimeDilationFactors()
        {
            if (_nodeTimeDilationFactor == null)
                InitAsynchronousTime();

            if (_correlationMass == null || _correlationMass.Length != N)
                RecomputeCorrelationMass();

            Parallel.For(0, N, i =>
            {
                // Compute local mass (Read-only access to _correlationMass and Weights)
                double localMass = _correlationMass[i];
                foreach (int j in Neighbors(i))
                {
                    if (j < _correlationMass.Length)
                        localMass += _correlationMass[j] * Weights[i, j];
                }

                // Compute curvature (Heavy calculation!)
                double localCurvature = 0;
                int edgeCount = 0;
                foreach (int j in Neighbors(i))
                {
                    localCurvature += CalculateGraphCurvature(i, j); // Ensure this method calculates purely and doesn't write
                    edgeCount++;
                }


                double avgCurvature = edgeCount > 0 ? localCurvature / edgeCount : 0;

                // Time dilation: dτ/dT = sqrt(max(ε, 1 - α*M - β*κ))
                const double alpha = 0.1;  // Mass coupling
                const double beta = 0.05;  // Curvature coupling
                const double epsilon = 0.01; // Minimum to prevent division by zero

                double factor = 1.0 - alpha * localMass - beta * avgCurvature;
                factor = Math.Max(epsilon, factor);
                factor = Math.Sqrt(factor);

                // Clamp to reasonable range
                _nodeTimeDilationFactor[i] = Math.Clamp(factor, MinTimeDilation, MaxTimeDilation);
            });
        }

        /// <summary>
        /// Perform asynchronous update step
        /// </summary>
        public void StepAsynchronous()
        {
            if (_nodeProperTime == null || _nodeNextUpdateTime == null || _nodeTimeDilationFactor == null)
                InitAsynchronousTime();

            // Update time dilation every step
            UpdateTimeDilationFactors();

            // Update nodes that are ready
            for (int i = 0; i < N; i++)
            {
                if (_nodeNextUpdateTime[i] <= _globalTime)
                {
                    // Compute proper time increment
                    double dTau = BaseTimestep * _nodeTimeDilationFactor[i];

                    // Update node state
                    UpdateSingleNode(i, dTau);

                    // Advance proper time and schedule next update
                    _nodeProperTime[i] += dTau;
                    _nodeNextUpdateTime[i] = _globalTime + BaseTimestep;
                }
            }

            // Advance global coordinate time
            _globalTime += BaseTimestep;
        }

        /// <summary>
        /// Update a single node's state
        /// </summary>
        private void UpdateSingleNode(int i, double dTau)
        {
            // Update node excitation state based on neighbors
            if (State[i] == NodeState.Refractory)
            {
                _refractoryCounter[i]--;
                if (_refractoryCounter[i] <= 0)
                {
                    State[i] = NodeState.Rest;
                }
                return;
            }

            // Count excited neighbors
            int excitedNeighbors = 0;
            double weightedExcitation = 0;

            foreach (int j in Neighbors(i))
            {
                if (State[j] == NodeState.Excited)
                {
                    excitedNeighbors++;
                    weightedExcitation += Weights[i, j];
                }
            }

            // Update based on local rules
            if (State[i] == NodeState.Excited)
            {
                // Transition to refractory
                State[i] = NodeState.Refractory;
                _refractoryCounter[i] = DynamicBaseRefractorySteps;

                // Propagate signals to neighbors (causal propagation)
                foreach (int j in Neighbors(i))
                {
                    PropagateSignal(i, j);
                }
            }
            else if (State[i] == NodeState.Rest)
            {
                // Activation probability based on neighbors
                double activationThreshold = 0.5;
                if (weightedExcitation > activationThreshold)
                {
                    State[i] = NodeState.Excited;
                }
            }

            // Update local fields if available
            if (LocalPotential != null && i < LocalPotential.Length)
            {
                // Diffusion and decay
                double diffusion = 0;
                int neighborCount = 0;

                foreach (int j in Neighbors(i))
                {
                    if (j < LocalPotential.Length)
                    {
                        diffusion += (LocalPotential[j] - LocalPotential[i]) * Weights[i, j];
                        neighborCount++;
                    }
                }

                if (neighborCount > 0)
                {
                    LocalPotential[i] += 0.1 * dTau * diffusion / neighborCount;
                }

                // Decay
                LocalPotential[i] *= (1.0 - 0.05 * dTau);
            }
        }


        public void UpdateHawkingProcess(double dt)
        {
            // 1. Цикл по всем ребрам (потенциальные места рождения пар)
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;

                    // 2. Вычисляем градиент гравитации (приливную силу)
                    // Разница кривизн на концах ребра
                    double R_i = GetLocalCurvature(i);
                    double R_j = GetLocalCurvature(j);
                    double tidalForce = Math.Abs(R_i - R_j);

                    // Порог рождения: чем выше приливная сила, тем вероятнее рождение пары
                    double productionProb = 0.01 * tidalForce * dt;

                    if (_rng.NextDouble() < productionProb)
                    {
                        // 3. Рождение пары (+E, -E)
                        // Частица с "минусом" падает туда, где кривизна (гравитация) больше
                        int heavyNode = (R_i > R_j) ? i : j;
                        int lightNode = (R_i > R_j) ? j : i;

                        // Энергия пары (из вакуума)
                        double pairEnergy = 0.1; // Mc^2

                        // "Излучение": улетает наружу (в область с меньшей кривизной)
                        if (LocalPotential != null) LocalPotential[lightNode] += pairEnergy;

                        // "Античастица": падает внутрь (уменьшает массу кластера)
                        if (_correlationMass != null) _correlationMass[heavyNode] -= pairEnergy;

                        // Баланс соблюден: +E + (-E) = 0.
                    }
                }
            }
        }

        /// <summary>
        /// Get node proper time
        /// </summary>
        public double GetNodeProperTime(int i)
        {
            if (_nodeProperTime == null || i < 0 || i >= N)
                return 0.0;

            return _nodeProperTime[i];
        }

        /// <summary>
        /// Get node time dilation factor
        /// </summary>
        public double GetNodeTimeDilation(int i)
        {
            if (_nodeTimeDilationFactor == null || i < 0 || i >= N)
                return 1.0;

            return _nodeTimeDilationFactor[i];
        }

        /// <summary>
        /// Get current global time
        /// </summary>
        public double GlobalTime => _globalTime;

        /// <summary>
        /// Compute local proper time step for a node based on its time dilation factor.
        /// Used by the event-based simulation to schedule next updates.
        /// RQ-Hypothesis Item 5: Local proper time τ varies with gravity.
        /// </summary>
        /// <param name="nodeId">Node index</param>
        /// <returns>Proper time increment dτ for this node</returns>
        public double ComputeLocalProperTimeStep(int nodeId)
        {
            if (nodeId < 0 || nodeId >= N)
                return BaseTimestep;

            // Time dilation from mass and curvature
            double timeDilation = GetNodeTimeDilation(nodeId);

            // Local proper time step: dτ = dt_base * time_dilation
            // Higher mass/curvature → smaller time dilation → slower local time
            double dTau = BaseTimestep * timeDilation;

            // Also consider local curvature for additional gravitational effect
            double localCurvature = GetLocalCurvature(nodeId);
            double curvatureFactor = 1.0 / Math.Sqrt(Math.Abs(localCurvature) + PhysicsConstants.CurvatureRegularizationEpsilon);
            curvatureFactor = Math.Clamp(curvatureFactor, 0.1, 10.0);

            dTau *= curvatureFactor;

            // Clamp to reasonable range
            return Math.Clamp(dTau, BaseTimestep * MinTimeDilation, BaseTimestep / MinTimeDilation);
        }

        /// <summary>
        /// Update physics for a single node and its neighbors with local proper time.
        /// RQ-Hypothesis Item 5: Each node evolves according to its own proper time τ.
        /// </summary>
        /// <param name="nodeId">Node to update</param>
        /// <param name="dTau">Local proper time step</param>
        public void UpdateSingleNodePhysics(int nodeId, double dTau)
        {
            if (nodeId < 0 || nodeId >= N)
                return;

            // Advance proper time for this node
            if (_nodeProperTime != null)
            {
                _nodeProperTime[nodeId] += dTau;
            }

            // === Update node state (excitation dynamics) ===
            UpdateSingleNode(nodeId, dTau);

            // === Update scalar field at this node ===
            if (ScalarField != null && _scalarMomentum != null)
            {
                // Local Klein-Gordon evolution
                double force = ComputeScalarFieldForce(nodeId);
                _scalarMomentum[nodeId] += dTau * force;
                ScalarField[nodeId] += dTau * _scalarMomentum[nodeId];
            }

            // === Update local potential ===
            if (LocalPotential != null && nodeId < LocalPotential.Length)
            {
                // Diffusion from neighbors (weighted by their time factors)
                double diffusion = 0.0;
                int neighborCount = 0;
                foreach (int j in Neighbors(nodeId))
                {
                    if (j < LocalPotential.Length)
                    {
                        double neighborTimeDilation = GetNodeTimeDilation(j);
                        diffusion += (LocalPotential[j] - LocalPotential[nodeId]) * Weights[nodeId, j] * neighborTimeDilation;
                        neighborCount++;
                    }
                }
                if (neighborCount > 0)
                {
                    LocalPotential[nodeId] += 0.1 * dTau * diffusion / neighborCount;
                }
            }

            // === Update edge weights connecting to this node (local geometry) ===
            foreach (int j in Neighbors(nodeId))
            {
                // Geometry update with proper time
                double curvature = CalculateGraphCurvature(nodeId, j);
                double deltaW = -dTau * PhysicsConstants.GravitationalCoupling * 0.01 * curvature;

                double newWeight = Weights[nodeId, j] + deltaW;
                newWeight = Math.Clamp(newWeight, 0.0, 1.0);
                Weights[nodeId, j] = newWeight;
                Weights[j, nodeId] = newWeight;
            }
        }

        /// <summary>
        /// Run a single step of event-based asynchronous simulation.
        /// Processes the next event from the priority queue and updates the node.
        /// RQ-Hypothesis Item 5: True relational time with local clocks.
        /// </summary>
        /// <returns>True if an event was processed, false if queue is empty</returns>
        public bool StepEventBased()
        {
            if (_eventQueue == null)
                InitAsynchronousTime();

            if (!_eventQueue!.TryDequeue(out var evt, out double time))
                return false;

            // Update global reference time
            CurrentReferenceTime = time;
            int nodeId = evt.NodeId;

            // Compute local proper time step for this node
            double dTau = ComputeLocalProperTimeStep(nodeId);

            // Update this node with its local proper time
            UpdateSingleNodePhysics(nodeId, dTau);

            // Schedule next update for this node
            double nextTime = time + dTau;
            _eventQueue.Enqueue(new NodeUpdateEvent(nextTime, nodeId, NodeEventType.Update), nextTime);

            return true;
        }

        /// <summary>
        /// Run event-based simulation for a batch of events.
        /// This is the main entry point for asynchronous time mode.
        /// RQ-Hypothesis Item 5: Priority queue-based local clock evolution.
        /// </summary>
        /// <param name="eventsPerStep">Number of events to process per call</param>
        public void StepEventBasedBatch(int eventsPerStep)
        {
            for (int i = 0; i < eventsPerStep; i++)
            {
                if (!StepEventBased())
                    break;
            }

            // Update time dilation factors periodically
            if (_eventQueue != null && _eventQueue.Count > 0)
            {
                UpdateTimeDilationFactors();
            }
        }
    }
}
