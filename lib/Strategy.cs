using System;
using System.Security.AccessControl;
namespace Phy.Lib;

public enum SIG
{
    HOLD = 101,   // 0
    BUY = 102,    // 1
    SELL = 103,   // 2
    CLOSE = 104, // 3
    TRADE = 105,// 4
    NOTRADE = 106  // 5
}


public interface IStrategy
{
    public SIG VolatilityMomentumSIG();
    public SIG GetSignal();

}
public class Strategy : IStrategy
{
    private PhysicsEngine _engine;

    private SIG _tacticalSignal;
    private Stats _stats;
    private Utils _utils;

    public Strategy(PhysicsEngine engine, Stats stats, Utils utils)
    {
        _engine = engine;
        _stats = stats;
        _utils = utils;
        _tacticalSignal = SIG.HOLD;
    }


    public SIG VolatilityMomentumSIG()
    {
        double strictness = 1.0; // Adjust this to make the strategy more or less aggressive
                                 // 1. KINETIC GATE (Is the market physically moving?)
                                 // Use your Universal Physics Engine. No more raw pip calculations.

        int SHIFT = _engine.GetIndData().Shift;

        double kineticEnergy = _engine.atrKinetic();

        // If Kinetic Energy is extremely low, the market is completely dead.
        if (kineticEnergy < (0.10 * strictness))
            return SIG.NOTRADE;

        // 2. POTENTIAL GATE (Is there a trend context?)
        // "Is the atmosphere charged?"
        double potential = _engine.adxPotential();

        // Gate: If potential is too low (< 0.75 = ADX 15), the market is wandering.
        if (potential < (0.75 * strictness))
            return SIG.NOTRADE;


        double normalizedScore = _engine.VolatilityEfficiency();

        // 4. DECISION: Dynamic Physics Threshold
        double baseThreshold = 0.40;

        // Discount Logic: Strong trends reduce the work required.
        // Lowered multiplier from 0.50 to 0.20 so we don't discount too aggressively.
        double discount = 0.20 * (potential - 1.0);
        double finalThreshold = (baseThreshold - discount) * strictness;

        // CRITICAL SAFETY FIX:
        // We MUST require at least *some* true volatility expansion (e.g., 0.15).
        // If we let this drop to 0.0 or negative, MathAbs() will allow false signals.
        finalThreshold = Math.Max(finalThreshold, 0.15);

        double slope30 = _stats.slopesVal(_engine.GetIndData().Ima30, shift: SHIFT)[SHIFT];

        // 5. TRIGGER
        if (Math.Abs(normalizedScore) > finalThreshold)
        {
            if (slope30 > 0.1) return SIG.BUY;
            else if (slope30 < -0.1) return SIG.SELL;
        }
        return SIG.NOTRADE;
    }
    public SIG GetSignal()
    {
        SIG signal = VolatilityMomentumSIG();
        _tacticalSignal = signal;
        return signal;
    }

    public SIG GetPhysicsSignal()
    {
        double velocity = _engine.GetVelocity(_engine.GetIndData().Ima30);
        double acceleration = _engine.GetAcceleration(_engine.GetIndData().Ima30);
        double zScore = _engine.GetMomentumZScore(_engine.GetIndData().Ima30);

        // If velocity is high, acceleration is positive, and it's not a 'freak' move (Z < 2.0)
        if (velocity > 0.5 && acceleration > 0 && zScore < 2.0)
        {
            return SIG.BUY;
        }

        // Mean Reversion: High velocity but negative acceleration at a high Z-score
        if (zScore > 2.5 && acceleration < 0)
        {
            return SIG.SELL;
        }

        return SIG.HOLD;
    }

}