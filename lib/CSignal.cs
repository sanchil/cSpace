using System;
using System.Security.AccessControl;
namespace Phy.Lib;

public interface ISignal
{
    public SIG GetPhysicsSignal();
    public SIG VolatilityMomentumSIG();
    public SIG GetSignal();
    public SIG TradeSlopeSIG(in DTYPE fast, in DTYPE slow, int magicnumber = -1);
    public SIG TradeSlopeSIG_Static(in DTYPE fast, in DTYPE slow, int magicnumber = -1);
    public SIG waveTideSIG(in DTYPE fast, in DTYPE med, in DTYPE slow);




}
public class CSignal : ISignal
{
    private PhysicsEngine _engine;

    private SIG _tacticalSignal;
    private CStats _stats;
    private CUtils _utils;

    private static readonly double[] closeRVal = { 1.3, 1.2, 1.1, 1.0, 0.9 };
    private double m_peakRatio;  // class member
    private DateTime m_last_bar;
    private SIG m_cached;


    public CSignal(PhysicsEngine engine, CStats stats, CUtils utils)
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

        double slope30 = _stats.slopesVal(_engine.GetIndData().Ima30, shift: SHIFT).val2;

        // 5. TRIGGER
        if (Math.Abs(normalizedScore) > finalThreshold)
        {
            if (slope30 > 0.1) return SIG.BUY;
            else if (slope30 < -0.1) return SIG.SELL;
        }
        return SIG.NOTRADE;
    }
    public T_SIG InitSignal()
    {
        T_SIG tSig = new T_SIG();
        IndData indData = _engine.GetIndData();
        tSig.volMomentumSIG = VolatilityMomentumSIG();
        tSig.tradeSlopeSIG = TradeSlopeSIG(_stats.slopesVal(indData.Ima30), _stats.slopesVal(indData.Ima60));
        tSig.slopeSIG = SlopeAnalyzerSIG(_stats.slopesVal(indData.Ima30));
        tSig.candleVolSIG = CandleVolSIG(indData.Open, indData.Close, indData.TickVolume, indData.Atr[indData.Shift]);
        tSig.singleCandleVolSIG = new SingleCandleVolSIG(_engine).Analyze(indData.Open, indData.Close, indData.TickVolume, indData.Atr[indData.Shift]);
        tSig.layeredMomentumSIG = LayeredMomentumSIG(indData.Ima30);
        tSig.physicsSIG = GetPhysicsSignal();

        return tSig;
    }

    //    public SIG GetSignal()
    //     {
    //         SIG signal = VolatilityMomentumSIG();
    //         _tacticalSignal = signal;
    //         return signal;
    //     }

    public SIG GetSignal()
    {
        T_SIG tSig = InitSignal();
        SIG signal = tSig.volMomentumSIG;
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


    public SIG TradeSlopeSIG(in DTYPE fast, in DTYPE slow, int magicnumber = -1)
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
        ReadOnlySpan<double> closeRatios = closeRVal;
        double CLOSERATIO = closeRatios[regimeIdx];

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

    public SIG TradeSlopeSIG_Static(in DTYPE fast, in DTYPE slow, int magicnumber = -1)
    {

        IndData indData = _engine.GetIndData();
        double atr = indData.Atr[indData.Shift];
        double floor = atr * 0.30;

        double fS = fast.val1;
        double sS = slow.val1;
        double absSlow = Math.Abs(sS);


        SIG dir = (sS > 0) ? SIG.BUY : SIG.SELL;
        //const double SLOPERANGELIMIT = 0.3;

        // 1. CALL THE PHYSICS ENGINE
        // This one call handles Directional Alignment AND the Structural Floor.
        //double mScore = ms.slopeRatio(fS, sS, SLOPERANGELIMIT);
        // To this:
        double mScore = _engine.expansionCompressionRatio(fS, sS, floor);


        // 2. THE LEAN POLICY
        // If mScore is 1.0, the Metric has already verified Direction, Floor, and Expansion.
        if (mScore >= 1.0) return dir;

        // For Case B (Compression), we check for the "Power Trend" extra requirement.
        if (mScore >= 0.8 && absSlow >= (floor * 1.5)) return dir;

        if (mScore == 0) return SIG.CLOSE;

        // If the Metric returned 0.0 (Veto) or a weak ratio, we bail.
        return SIG.CLOSE;
    }

    public SIG microWaveSIG(in DTYPE fast, in DTYPE med)
    {

        // 1. THE MICRO-FLOOR (Aggressive)
        // We lower the barrier to entry. We only need 10% of ATR to consider it "Active."
        IndData indData = _engine.GetIndData();
        double atr = indData.Atr[indData.Shift];

        double MICRO_FLOOR = atr * 0.10;

        double fS = fast.val1;
        double mS = med.val1;
        SIG dir = (fS > 0) ? SIG.BUY : SIG.SELL;

        // 2. THE VELOCITY CHECK (The "Explosion" Gate)
        // We use slopeRatio but we pass our lower MICRO_FLOOR.
        // We want to see the Wave pulling away from the Current.
        double vScore = _engine.expansionCompressionRatio(fS, mS, MICRO_FLOOR);

        // 3. THE MICRO-POLICY
        // We ONLY enter if the expansion is nearly perfect (>= 0.95)
        // This ensures we are catching the "Meat" of the micro-move.
        if (vScore >= 0.95)
        {
            return dir;
        }

        // 4. THE LIGHTNING EXIT
        // If the velocity score drops even slightly (e.g., below 0.70),
        // we BAIL. There is no macro structure to save us here.
        if (vScore < 0.70) return SIG.CLOSE;

        return SIG.NOSIG;
    }

    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public SIG waveTideSIG(in DTYPE fast, in DTYPE med, in DTYPE slow)
    {

        // THE TRIPLE-GEOMETRY CHAIN
        SIG waveSignal = TradeSlopeSIG_Static(fast, med);  // Micro-Expansion
        SIG tideSignal = TradeSlopeSIG_Static(med, slow);  // Macro-Expansion

        if (waveSignal == SIG.BUY && tideSignal == SIG.BUY) return SIG.BUY;
        if (waveSignal == SIG.SELL && tideSignal == SIG.SELL) return SIG.SELL;
        return SIG.NOSIG;
    }

    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public SIG SlopeAnalyzerSIG(in DTYPE slope)
    {
        SlopeAnalyzer analyzer = new SlopeAnalyzer(_engine);
        return analyzer.Analyze(slope.val1);
    }

    //+------------------------------------------------------------------+
    //| Layered Filter: ADX → Histogram for Momentum Strength            |
    //+------------------------------------------------------------------+
    SIG LayeredMomentumSIG(in double[] signal, int N = 20)
    {

        double gate = _engine.layeredMomentumFilter(signal, N);
        if (gate == 0)
            return SIG.NOSIG;
        if (gate == 1)
            return SIG.BUY;
        if (gate == -1)
            return SIG.SELL;
        return SIG.NOSIG;
    }


    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    SIG CandleVolSIG(
            in double[] open,
            in double[] close,
            in double[] volume,
            in double atr,
            int period = 30,
            int SHIFT = 1
            )
    {

        IndData indData = _engine.GetIndData();
        double atr_pips = atr / indData.Point;

        double atr_floor = atr_pips * 0.05;

        // 1. DYNAMIC WINDOWS (Guaranteed to be distinct)
        //int fast_n = (int)(period * 0.5);
        int fast_n = Math.Max(5, (int)(period * 0.4));  // 40% instead of 50% — faster response

        if (fast_n < 5) fast_n = 5; // Minimum 5 candles to avoid noise

        // 2. RAW KINEMATICS
        double slow = _engine.vWCM_Raw(10);
        double fast = _engine.vWCM_Raw(10);

        if (Math.Abs(slow) < atr_floor) return SIG.NOSIG;  // veto very flat volume


        // 3. DIRECTIONAL AGREEMENT
        bool agree_dir = (slow > 0 && fast > 0) || (slow < 0 && fast < 0);

        // 4. MOMENTUM RATIO (Fast must maintain at least 75% of Slow's power)
        // We use MathAbs to safely compare the force regardless of direction
        //bool agree_str = (MathAbs(fast) / (MathAbs(slow) + DBL_EPSILON)) > 0.75;
        bool agree_str = (Math.Abs(fast) / (Math.Abs(slow) + indData.DBL_EPSILON)) > 0.6;

        // 5. THE VERDICT
        SIG sig = SIG.NOSIG;
        if (agree_dir && agree_str)
        {
            sig = (slow > 0) ? SIG.BUY : SIG.SELL;
        }

        // [Optional Debugging Log]
        // PrintFormat("vWCM | Slow:%.4f Fast:%.4f → %s", slow, fast, util.getSigString(sig));

        return sig;
    }

}



//+------------------------------------------------------------------+
//| singleCandleVolSIG - Final, bulletproof version                     |
//+------------------------------------------------------------------+

class SingleCandleVolSIG
{
    private readonly IPhysicsEngine _engine;

    private DateTime last_bar;
    private SIG cached;

    public SingleCandleVolSIG(IPhysicsEngine eng)
    {
        _engine = eng;
        last_bar = DateTime.MinValue; ;
        cached = SIG.NOSIG;
    }
    public SIG Analyze(
       in double[] open,
       in double[] close,
       in double[] volume,
       in double atr,
       int period = 30,
       int SHIFT = 1
       )
    {
        IndData indData = _engine.GetIndData();


        if (indData.Time[0] == last_bar)
            return cached;


        last_bar = indData.Time[0];

        double atr_pips = atr / indData.Point;
        //if(atr_pips < 8.0) {
        //   cached = SAN_SIGNAL::NOSIG;
        //   return cached;
        //}
        //double slow = stats.vWCM_Score(open, close, volume, period,0,SHIFT);
        double slow = _engine.vWCM_Raw(10);
        // Print("[SLOWVCM]: " + slow);

        if ((slow > -0.05) && (slow < 0.1))
            cached = SIG.NOSIG;
        if (slow >= 0.1)
            cached = SIG.BUY;
        if (slow <= -0.05)
            cached = SIG.SELL;

        //cached = (slow > 0) ? SAN_SIGNAL::BUY : SAN_SIGNAL::SELL;

        //PrintFormat("vWCM | ATR:%.1f pips | Slow:%.4f",
        //            atr_pips, slow,
        //            cached==BUY?"BUY":cached==SELL?"SELL":"NOSIG");

        return cached;
    }


}





class SlopeAnalyzer
{
    // --- 1. State Memory (Moved from Function-Static to Class-Fields) ---
    // These replace 'static double peakPositive = 0;'
    private double _peakPositive = 0;
    private double _peakNegative = 0;
    private SIG _currentIdx = SIG.NOTRADE;

    private readonly IPhysicsEngine _engine;

    public SlopeAnalyzer(IPhysicsEngine engine)
    {
        _engine = engine;
    }

    public SIG Analyze(double slopeValue)
    {
        // Constants
        const double BASE_DECAY = 0.8;
        const double MIN_SLOPE = 0.2;
        const double HYSTERESIS = 0.90;

        double s = slopeValue;

        // 2. Adaptive Logic (Using your ported Physics Engine)
        double adxNorm = _engine.adxKinetic();
        double adaptedDecay = BASE_DECAY + (0.18 * adxNorm);

        // --- LOGIC GATE ---

        // A. RESET LOGIC
        if (s < -MIN_SLOPE) _peakPositive = 0;
        if (s > MIN_SLOPE) _peakNegative = 0;

        // B. BUY LOGIC
        double buyThreshold = (_peakPositive > 0) ? (_peakPositive * adaptedDecay) : MIN_SLOPE;

        if (s > buyThreshold)
        {
            _peakPositive = Math.Max(_peakPositive, s);
            _currentIdx = SIG.BUY;
            return SIG.BUY;
        }

        // C. SELL LOGIC
        double sellThreshold = (_peakNegative < 0) ? (_peakNegative * adaptedDecay) : -MIN_SLOPE;

        if (s < sellThreshold)
        {
            _peakNegative = Math.Min(_peakNegative, s);
            _currentIdx = SIG.SELL;
            return SIG.SELL;
        }

        // D. EXIT LOGIC with HYSTERESIS
        if (_currentIdx == SIG.BUY)
        {
            double exitLevel = (_peakPositive * adaptedDecay) * HYSTERESIS;
            if (s < exitLevel)
            {
                _currentIdx = SIG.NOTRADE;
                return SIG.CLOSE;
            }
            return SIG.BUY; // HOLD
        }

        if (_currentIdx == SIG.SELL)
        {
            double exitLevel = (_peakNegative * adaptedDecay) * HYSTERESIS;
            if (s > exitLevel)
            {
                _currentIdx = SIG.NOTRADE;
                return SIG.CLOSE;
            }
            return SIG.SELL; // HOLD
        }

        return SIG.NOTRADE;
    }
}