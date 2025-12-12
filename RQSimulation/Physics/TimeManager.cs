using System;

namespace RQSimulation.Physics
{
    /// <summary>
    /// Provides unified relational time stepping. Computes the time
    /// increment dTau from the system's internal clock subsystem and
    /// exposes a single value to be used for all physics updates. This
    /// resolves checklist item #3 by removing global fixed dt.
    /// </summary>
    public static class TimeManager
    {
        /// <summary>
        /// Computes the relational time step dTau based on a set of
        /// clock amplitudes. The result is the norm of the difference
        /// between successive clock states. The returned value is
        /// clamped to a reasonable range [min, max] to avoid unstable
        /// integration.
        /// </summary>
        /// <param name="previousClockState">State of the clock at the previous step.</param>
        /// <param name="currentClockState">State of the clock at the current step.</param>
        /// <param name="minStep">Minimum allowed time step.</param>
        /// <param name="maxStep">Maximum allowed time step.</param>
        /// <returns>The relational time increment dTau.</returns>
        public static double ComputeRelationalTimeStep(
            double[] previousClockState,
            double[] currentClockState,
            double minStep = 1e-3,
            double maxStep = 1e-1)
        {
            ArgumentNullException.ThrowIfNull(previousClockState);
            ArgumentNullException.ThrowIfNull(currentClockState);

            if (previousClockState.Length != currentClockState.Length)
            {
                throw new ArgumentException("Clock states must be of equal length.");
            }

            // Euclidean norm of the difference between clock states
            double sumSq = 0.0;
            for (int i = 0; i < previousClockState.Length; i++)
            {
                double delta = currentClockState[i] - previousClockState[i];
                sumSq += delta * delta;
            }
            double norm = Math.Sqrt(sumSq);

            // Clamp to [minStep, maxStep] to avoid extremes
            if (double.IsNaN(norm) || double.IsInfinity(norm))
                norm = minStep;

            return Math.Clamp(norm, minStep, maxStep);
        }
    }
}
