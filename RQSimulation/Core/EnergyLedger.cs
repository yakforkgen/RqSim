using System;

namespace RQSimulation
{
    /// <summary>
    /// Tracks energy conservation across all simulation operations
    /// Ensures strict adherence to energy conservation laws
    /// Implements RQ-hypothesis checklist item 4: Unified Hamiltonian and Energy Ledger
    /// 
    /// Energy accounting:
    /// - VacuumPool: Available vacuum energy for topology changes
    /// - MatterEnergy: Energy in particle clusters
    /// - FieldEnergy: Energy in gauge fields
    /// - Total = VacuumPool + MatterEnergy + FieldEnergy (conserved)
    /// </summary>
    public class EnergyLedger
    {
        private double _totalEnergy;
        private double _externalInjection;
        private double _vacuumBorrowing;
        private double _vacuumPool;
        private double _matterEnergy;
        private double _fieldEnergy;
        private const double Tolerance = 1e-6;
        private bool _initialized;

        /// <summary>
        /// Vacuum energy pool available for topology changes and particle creation.
        /// Implements checklist item 4.1: Unified energy functional.
        /// </summary>
        public double VacuumPool 
        { 
            get => _vacuumPool;
            private set => _vacuumPool = value;
        }
        
        /// <summary>
        /// Energy stored in matter clusters
        /// </summary>
        public double MatterEnergy => _matterEnergy;
        
        /// <summary>
        /// Energy stored in gauge fields
        /// </summary>
        public double FieldEnergy => _fieldEnergy;

        /// <summary>
        /// Initialize the ledger with the current total energy
        /// </summary>
        public void Initialize(double initialEnergy)
        {
            _totalEnergy = initialEnergy;
            _externalInjection = 0;
            _vacuumBorrowing = 0;
            _vacuumPool = initialEnergy * PhysicsConstants.InitialVacuumPoolFraction;
            _matterEnergy = 0;
            _fieldEnergy = 0;
            _initialized = true;
        }
        
        /// <summary>
        /// Initialize with detailed energy breakdown
        /// </summary>
        public void Initialize(double vacuumEnergy, double matterEnergy, double fieldEnergy)
        {
            _vacuumPool = vacuumEnergy;
            _matterEnergy = matterEnergy;
            _fieldEnergy = fieldEnergy;
            _totalEnergy = vacuumEnergy + matterEnergy + fieldEnergy;
            _externalInjection = 0;
            _vacuumBorrowing = 0;
            _initialized = true;
        }

        /// <summary>
        /// Try to spend energy from vacuum pool for topology change.
        /// Returns false if not enough energy is available.
        /// Implements checklist item 4.3: Energy check for topology changes.
        /// </summary>
        /// <param name="amount">Energy required for the change</param>
        /// <returns>True if energy was available and spent, false otherwise</returns>
        public bool TrySpendVacuumEnergy(double amount)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized. Call Initialize() first.");
            }
            
            if (amount < 0)
            {
                throw new ArgumentException("Cannot spend negative energy", nameof(amount));
            }
            
            if (_vacuumPool >= amount)
            {
                _vacuumPool -= amount;
                return true;
            }
            return false; // Insufficient energy - action forbidden
        }
        
        /// <summary>
        /// Return energy to vacuum pool (e.g., from radiation, cluster decay).
        /// Implements checklist item 4.3: Energy return mechanism.
        /// </summary>
        /// <param name="amount">Energy to return to vacuum</param>
        public void RegisterRadiation(double amount)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized. Call Initialize() first.");
            }
            
            if (amount < 0)
            {
                throw new ArgumentException("Cannot register negative radiation", nameof(amount));
            }
            
            _vacuumPool += amount;
        }
        
        /// <summary>
        /// Check if a topology change can be afforded energy-wise.
        /// Implements checklist item 4.3: Metropolis criterion energy check.
        /// </summary>
        /// <param name="deltaEnergy">Energy change (positive = costs energy)</param>
        /// <returns>True if change is allowed, false otherwise</returns>
        public bool CanAfford(double deltaEnergy)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized. Call Initialize() first.");
            }
            
            if (deltaEnergy <= 0)
                return true; // Energy-releasing changes are always allowed
                
            return _vacuumPool >= deltaEnergy;
        }
        
        /// <summary>
        /// Update matter energy tracking
        /// </summary>
        public void UpdateMatterEnergy(double newMatterEnergy)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized. Call Initialize() first.");
            }
            
            double delta = newMatterEnergy - _matterEnergy;
            _matterEnergy = newMatterEnergy;
            
            // Compensate vacuum pool to maintain total energy conservation
            _vacuumPool -= delta;
        }
        
        /// <summary>
        /// Update field energy tracking
        /// </summary>
        public void UpdateFieldEnergy(double newFieldEnergy)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized. Call Initialize() first.");
            }
            
            double delta = newFieldEnergy - _fieldEnergy;
            _fieldEnergy = newFieldEnergy;
            
            // Compensate vacuum pool to maintain total energy conservation
            _vacuumPool -= delta;
        }
        
        /// <summary>
        /// Get total tracked energy (should be constant)
        /// </summary>
        public double TotalTrackedEnergy => _vacuumPool + _matterEnergy + _fieldEnergy;

        /// <summary>
        /// Record energy injection from external source (e.g., impulse)
        /// </summary>
        public void RecordExternalInjection(double energy, string source)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized. Call Initialize() first.");
            }

            _externalInjection += energy;
            _vacuumPool += energy; // External energy goes to vacuum pool
            Console.WriteLine($"[ENERGY] External injection: {energy:F6} from {source}");
        }

        /// <summary>
        /// Record energy borrowed from vacuum (must be repaid)
        /// </summary>
        public void BorrowFromVacuum(double energy)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized.");
            }

            _vacuumBorrowing += energy;
        }

        /// <summary>
        /// Repay energy borrowed from vacuum
        /// </summary>
        public void RepayToVacuum(double energy)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized.");
            }

            _vacuumBorrowing -= energy;
            if (_vacuumBorrowing < -1e-10) // Allow small numerical error
            {
                throw new InvalidOperationException(
                    $"Vacuum debt cannot be negative: {_vacuumBorrowing:F10}. " +
                    "More energy repaid than borrowed.");
            }
        }

        /// <summary>
        /// Validate that energy is conserved within tolerance
        /// </summary>
        public void ValidateConservation(double currentEnergy)
        {
            if (!_initialized)
            {
                throw new InvalidOperationException("EnergyLedger not initialized.");
            }

            double expected = _totalEnergy + _externalInjection;
            double error = Math.Abs(currentEnergy - expected);
            double relativeError = Math.Abs(expected) > 1e-10 
                ? error / Math.Abs(expected) 
                : error;

            if (error > Tolerance && relativeError > Tolerance)
            {
                throw new EnergyConservationException(
                    $"Energy conservation violated!\n" +
                    $"Expected: {expected:F8}\n" +
                    $"Current:  {currentEnergy:F8}\n" +
                    $"Error:    {error:F8} ({relativeError * 100:F2}%)\n" +
                    $"Injected: {_externalInjection:F8}\n" +
                    $"Vacuum:   {_vacuumBorrowing:F8}");
            }

            // Update total energy to current value (for next validation)
            _totalEnergy = currentEnergy;
        }

        /// <summary>
        /// Get current vacuum borrowing (should be ~0 in steady state)
        /// </summary>
        public double VacuumDebt => _vacuumBorrowing;

        /// <summary>
        /// Get total external energy injected
        /// </summary>
        public double TotalExternalInjection => _externalInjection;

        /// <summary>
        /// Get current tracked total energy
        /// </summary>
        public double TrackedEnergy => _totalEnergy;

        /// <summary>
        /// Reset external injection counter (e.g., at start of run)
        /// </summary>
        public void ResetExternalInjection()
        {
            _externalInjection = 0;
        }
    }

    /// <summary>
    /// Exception thrown when energy conservation is violated
    /// </summary>
    public class EnergyConservationException : Exception
    {
        public EnergyConservationException(string message) : base(message) { }
    }
}
