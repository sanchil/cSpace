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
    NOTRADE = 106,
    NOSIG = 107 // 5

}

public struct DTYPE
{
    public double val1;
    public double val2;
    public double val3;
    public double val4;
    public double val5;
};

public interface ISignal
{
    public SIG GetPhysicsSignal();
    public SIG VolatilityMomentumSIG();
    public SIG GetSignal();
    public SIG TradeSlopeSIG(in DTYPE fast, in DTYPE slow,int magicnumber = -1);
    
    


}
public class CSignal : ISignal
{
    private PhysicsEngine _engine;

    private SIG _tacticalSignal;
    private Stats _stats;
    private Utils _utils;

    private static readonly double[] closeRVal = { 1.3, 1.2, 1.1, 1.0, 0.9 };
    private double m_peakRatio;  // class member
    private DateTime m_last_bar;
    private SIG m_cached;


    public CSignal(PhysicsEngine engine, Stats stats, Utils utils)
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


    public SIG TradeSlopeSIG(in DTYPE fast, in DTYPE slow,int magicnumber = -1)
    {

        // --- 1. BAR OPENING CHECK (Fixed State Management) ---
        // Uses class members m_last_bar and m_cached instead of static variables
        // Passed Time[0] as a parameter or fetched cleanly.
        IndData indData = _engine.GetIndData();
        if (indData.Time[0] == m_last_bar)
            return m_cached;

        m_last_bar = indData.Time[0];
        double atr = indData.Atr[indData.Shift];

        // --- 2. DYNAMIC CONSTANTS & INPUTS ---
        const double MIN_SLOW_THRESHOLD = 0.0001; // Avoid division by zero

        // The Elastic Ruler: We define slope thresholds as fractions of the ATR.
        // Example: If ATR is 10 pips, FLAT is anything moving less than 0.5 pips/bar.
        double FLAT_REGIME_ENTRY = atr * 0.05;  // 5% of ATR
        double WEAK_TREND_BOUND = atr * 0.10;  // 10% of ATR
        double MID_TREND_BOUND = atr * 0.20;  // 20% of ATR
        double STRONG_TREND_BOUND = atr * 0.35;  // 35% of ATR
        double TURBO_TREND_BOUND = atr * 0.50;  // 50% of ATR

        double fastSlope = fast.val1;
        double slowSlope = slow.val1;
        double absSlow = Math.Abs(slowSlope);

        // --- 3. ADAPTIVE REGIME SELECTION ---
        // We replace your static 0.35, 0.80, 1.50 guesses with the ATR boundaries.
        int regimeIdx = (absSlow <= WEAK_TREND_BOUND) ? 0 :
                        (absSlow <= MID_TREND_BOUND) ? 1 :
                        (absSlow <= STRONG_TREND_BOUND) ? 2 :
                        (absSlow <= TURBO_TREND_BOUND) ? 3 : 4;

        // .NET 8 optimized: The compiler handles the safety and performance automatically
        // ReadOnlySpan<double> closeRVal = [1.3, 1.2, 1.1, 1.0, 0.9];
        // This lives in memory once and never moves. Very fast.
        double CLOSERATIO = closeRVal[regimeIdx];

        double PEAK_DROP = _engine.getVolAdaptiveRetention();
        PEAK_DROP = Math.Max(Math.Min(PEAK_DROP, 0.99), 0.70); // Hard clamp

        // --- 4. CORE LOGIC BRANCHING ---

        // BRANCH A: The "Flat" Market (Singularity Avoidance)
        if (absSlow < MIN_SLOW_THRESHOLD)
        {
            // Reversal Check
            if ((m_cached == SIG.BUY && fastSlope < -FLAT_REGIME_ENTRY) ||
                  (m_cached == SIG.SELL && fastSlope > FLAT_REGIME_ENTRY))
            {
                m_peakRatio = 0;
                m_cached = SIG.CLOSE;
                return SIG.CLOSE;
            }

            // Entry Logic
            if (Math.Abs(fastSlope) > FLAT_REGIME_ENTRY)
            {
                if (m_peakRatio == 0) m_peakRatio = CLOSERATIO * 1.05;
                m_cached = (fastSlope > 0) ? SIG.BUY : SIG.SELL;
                return m_cached;
            }

            // FLAT MARKET ORPHAN FIX: If momentum died while holding a trade, close it.
            if (m_peakRatio > 0)
            {
                m_peakRatio = 0;
                m_cached = SIG.CLOSE;
                return SIG.CLOSE;
            }

            m_peakRatio = 0;
            m_cached = SIG.NOSIG;
            return SIG.NOSIG;
        }

        // BRANCH B: The Standard Adaptive Engine
        double ratio = fastSlope / slowSlope;

        // // Use PrintFormat for efficiency
        // PrintFormat("Ratio=%.3f | Peak=%.3f | DropLimit=%.3f | Regime=%d",
        //             ratio, m_peakRatio, (PEAK_DROP * m_peakRatio), regimeIdx);

        // 1. INSTANT REVERSAL (Divergence Check)
        if (ratio <= 0)
        {
            m_peakRatio = 0;
            m_cached = SIG.CLOSE;
            return SIG.CLOSE;
        }

        // 2. MOMENTUM DECAY EXIT (The Adaptive Stop)
        if (m_peakRatio > 0 && ratio < (PEAK_DROP * m_peakRatio))
        {
            m_peakRatio = 0;
            m_cached = SIG.CLOSE;
            return SIG.CLOSE;
        }

        // 3. WEAK ALIGNMENT EXIT (Hard Floor)
        if (ratio <= CLOSERATIO)
        {
            m_peakRatio = 0;
            m_cached = SIG.CLOSE;
            return SIG.CLOSE;
        }

        // 4. ENTRY & CONTINUATION (New Peak Tracking)
        if (ratio > CLOSERATIO)
        {
            if (ratio > m_peakRatio) m_peakRatio = ratio;
            m_cached = (fastSlope > 0) ? SIG.BUY : SIG.SELL;
            return m_cached;
        }

        m_cached = SIG.NOSIG;
        return SIG.NOSIG;
    }

}