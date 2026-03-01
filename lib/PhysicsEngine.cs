using System;
using System.Net.NetworkInformation;
using MathNet.Numerics.Statistics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace Phy.Lib;

public interface IBotEngine
{
    public void InitIndData();
    public double[] GetHistory(Func<int, double> getPricePtr, int count, int shift = 1);
    public void printData(in IndData data);
}

public interface IPhysicsEngine
{
    public double atrKinetic();
    public double adxKinetic(double scale = 50.0, int shift = 1);
    public double adxPotential(int period = 14);
    public double adxVector();
    public double VolatilityEfficiency();
    public double efficiencyRatio(in double[] sig, int period = 14);
    public double vWCM_Raw(int N = 10);
    public double vWCM_Smooth(int N = 10);
    public double volatilityAnomaly();
    public bool isTrendAccelerating(in double[] sig, int shift = 1);
    public double trendAccelStrength(in double[] sig, int shift = 1);
    public double trendQuality(in double[] sig, int period = 14);

    public double expansionCompressionRatio(in double fastS, in double slowS, in double SLOPEFLOOR = 0.3);
    public double springForce(double currentPrice, double smaValue);
    public double slopeAccelerationRatio(in double fSlope, in double mSlope, in double sSlope);
    public double geometricFanScore(in double fastS, in double medS, in double slowS);


    public double fractalAlignment(in double fastS, in double medS, in double slowS);
    public double kinematicAcceleration(in double fastS, in double slowS, in double slopeFloor = 0.3);

    public double overallMarketForce(int period);
    public double layeredMomentumFilter(in double[] values, int N = 20);

    public double GetVelocity(double[] sig, int period = 10);
    public double GetAcceleration(double[] sig, int period = 10);
    public double GetMomentumZScore(double[] sig, int period = 20);

    public double getLinearTimeRetention(int barsHeld, double decayRate = 0.05, double floor = 0.60);

    public double getVolAdaptiveRetention();
    public double getHybridRetention(int barsHeld);
    public double getHybridRetention_v2(int barsHeld, double trendQualityScore = 0.0);
    public double getPeakDecay(DECAY_STRATEGY strat = DECAY_STRATEGY.STRAT_ATR, double period = 14, int shift = 1);

    public double bayesianNextBarProb(bool volHigh, bool accelStrong);
    public FEATURE_VECTOR getFeatureVector(in IndData indData, in int SHIFT = 1);

    public double marketIntensity(in FEATURE_VECTOR fV);
    public double marketRegime(in FEATURE_VECTOR fV);
    public double neuronHoldScore(
       in double[] maArray, in double[] closeArray, in double[] openArray, in double[] volumeArray,
       int barsHeld, double atr, bool useOverallForce = true);
    public double bayesianHoldScore(
   in double[] maArray,
   in double[] closeArray,
   in double[] openArray,
   in double[] volumeArray,
   int barsHeld,
   double atr,
   bool useOverallForce = true
   );

    public int getHyperbolicCombinedScore(double b, double n, double f, double fra);
    public int getCobbDouglasCombinedScore(double b, double n, double f, double fra);



}

public readonly record struct IndData
{
    // --- 1. HISTORICAL ARRAYS (The "Primary" 500 & "Secondary" 120) ---
    // Note: In C#, arrays inside a struct are Reference Types pointing to the Heap.
    public double[] Open { get; init; }
    public double[] High { get; init; }
    public double[] Low { get; init; }
    public double[] Close { get; init; }
    public DateTime[] Time { get; init; }
    public double[] TickVolume { get; init; }
    // Indicators
    public double[] StdClose { get; init; }
    public double[] StdOpen { get; init; }
    public double[] Mfi { get; init; }
    public double[] Obv { get; init; }
    public double[] Rsi { get; init; }
    public double[] Atr { get; init; }
    public double[] Adx { get; init; }
    public double[] AdxPlus { get; init; }
    public double[] AdxMinus { get; init; }
    public double[] Ima5 { get; init; }
    public double[] Ima14 { get; init; }
    public double[] Ima30 { get; init; }
    public double[] Ima60 { get; init; }
    public double[] Ima120 { get; init; }
    public double[] Ima240 { get; init; }
    public double[] Ima500 { get; init; }
    public double[] AvgStd { get; init; }

    // --- 2. SCALARS (Trading State) ---
    public ulong MagicNumber { get; init; }
    public double CloseProfit { get; init; }
    public double StopLoss { get; init; }
    public double CurrProfit { get; init; }
    public double MaxProfit { get; init; }
    public SIG TradePosition { get; init; }
    public int CurrSpread { get; init; }
    public int Shift { get; init; }
    public int BarsHeld { get; init; }

    // --- 3. PHYSICS SCORES ---
    public double HoldScore { get; init; }
    public double BayesianHoldScore { get; init; }
    public double NeuronHoldScore { get; init; }
    public double FMSR_Norm { get; init; }
    public double FractalAlignment { get; init; }
    public double MicroLots { get; init; }
    public double ConvictionFactor { get; init; }
    public int PhysicsAction { get; init; }
    public double PipValue { get; init; }

    public double Point { get; init; }
    public long _Period { get; init; }
    public double PipSize { get; init; }

    public double Current_Period { get; init; }
    public double DBL_EPSILON { get; init; }
    public int Digits { get; init; }
}

public struct FEATURE_VECTOR
{

    // Group A: Momentum & Velocity (The Drivers)
    public double slopeIma5 { get; set; }
    public double slopeIma30 { get; set; }
    public double adxPlusMinusDiff { get; set; }
    public double rsi { get; set; }

    // Group B: Volatility & Energy (The Fuel)
    public double atr { get; set; }
    public double stdDevCP { get; set; }
    public double adx { get; set; }
    public double tVol { get; set; }

    // Group C: Structure & Stretch (The Geometry)

    public double priceElasticity { get; set; }
    public double mfi { get; set; }
    public double vWCM { get; set; }
    public double expansionCompression { get; set; }

    // Group D: The "Super-Feature" Additions.
    public double bayesianScore { get; set; }
    public double neuronScore { get; set; }
    public double fMSR_Norm { get; set; }
    public double fractalAlignment { get; set; }

};

public enum DECAY_STRATEGY
{
    STRAT_ATR,       // Volatility fuel (your current)
    STRAT_ADX,       // Trend quality
    STRAT_ER,        // Efficiency
    STRAT_MIX        // Weighted mix (0.5 ATR + 0.3 ADX + 0.2 ER)
};

public class PhysicsEngine : IPhysicsEngine
{
    private IndData _indData;
    private Stats _stats;
    private Utils _utils;

    private int SHIFT = 0;
    // Constructor



    public PhysicsEngine(IndData indData, Stats stats, Utils utils)
    {
        _indData = indData;
        _stats = stats;
        _utils = utils;
        SHIFT = _indData.Shift;
    }

    // Method-based initialization (Setter)
    public void SetIndData(IndData data)
    {
        _indData = data;
    }

    public IndData GetIndData() => _indData;

    // 4. atrKinetic (Universal Timeframe Logic - Sqrt Rule)
    // TRUTH: "Is the absolute movement large enough to trade?"
    // RETURNS: 0.0 to 1.0 (Normalized by Square Root of Time)
    public double atrKinetic()
    {
        double pipValue = _indData.PipSize;
        double _Period = _indData.Current_Period;
        double atr = _indData.Atr[SHIFT];


        // Safety check for weird broker data
        if (pipValue <= 0)
            return 0.0;

        double atrPips = atr / pipValue;

        // 1. BASELINE CALIBRATION (M15)
        // This is the "Energy Reference".
        // 30 pips on M15 is considered "Max Energy" (Score 1.0).
        // If you want to require MORE movement to trigger, increase this to 35.0.
        double baseRef = 30.0;

        // 2. PHYSICS SCALING: Square Root Rule
        // Formula: Ref_Current = Ref_Base * Sqrt(Period / BasePeriod)
        // We anchor to M15 (15 minutes).
        // M15 Ratio = 1.0. H1 Ratio = 2.0 (Sqrt 4). H4 Ratio = 4.0 (Sqrt 16).
        double timeRatio = (double)_Period / 15.0;

        // Prevent div/0 or negative sqrt
        if (timeRatio <= 0)
            timeRatio = 1.0;

        double physicsCeiling = baseRef * Math.Sqrt(timeRatio);

        // 3. NORMALIZE
        double atrNorm = Math.Min(Math.Max(atrPips / physicsCeiling, 0.0), 1.0);

        // 4. SQUASH (Kinetic Energy = v^2)
        // Punish weak moves, reward strong ones.
        return (atrNorm * atrNorm);
    }

    // 1. adxPotential (The Fuel Gauge)
    // TRUTH: "How charged is the market environment?"
    // LOGIC: Hybrid. Exponential below 20 (suppress noise), Linear above 20 (preserve scale).
    // RETURNS: 0.0 to 3.0+ (1.0 = Baseline ADX 20)
    public double adxPotential(int period = 14)
    {
        double adx = _indData.Adx[SHIFT];
        double raw = adx / 20.0;
        // Smart Hybrid Curve: Square the noise, keep the trend linear.
        return (raw < 1.0) ? (raw * raw) : raw;
    }

    public double GetVelocity(double[] sig, int period = 10)
    {
        // Velocity is the average slope over the last N bars
        // We use the slopesVal we already built
        double[] slopes = _stats.slopesVal(sig, SLOPEDENOM: period, shift: _indData.Shift);
        return slopes[0]; // Primary slope
    }

    public double GetAcceleration(double[] sig, int period = 10)
    {
        int S = _indData.Shift;
        // Acceleration = (Current Velocity - Previous Velocity) / Time
        double v1 = GetVelocity(sig, period); // Current 

        // To get previous velocity, we shift the calculation by 1
        double[] prevSlopes = _stats.slopesVal(sig, SLOPEDENOM: period, shift: S + 1);
        double v0 = prevSlopes[0];

        return v1 - v0;
    }

    public double GetMomentumZScore(double[] sig, int period = 20)
    {
        // Here we use MathNet via our Stats class
        // We want to know if the current 'Velocity' is an outlier
        double[] velocityHistory = new double[period];

        for (int i = 0; i < period; i++)
        {
            velocityHistory[i] = GetVelocity(sig, i + _indData.Shift);
        }

        var distribution = _stats.GetDistribution(velocityHistory, 0);
        return distribution.zScore;
    }

    public double VolatilityEfficiency()
    {
        double stdCp = _indData.StdClose[SHIFT];
        double stdOpen = _indData.StdOpen[SHIFT];
        double slopeCP = _stats.slopesVal(_indData.StdClose, shift: _indData.Shift)[SHIFT];
        double slopeOP = _stats.slopesVal(_indData.StdOpen, shift: _indData.Shift)[SHIFT];

        // A. Structure Ratio (Current / Open)
        double denominator = (stdOpen < 0.00005) ? 0.00005 : stdOpen;
        double structureRatio = stdCp / denominator;

        // B. Momentum Ratio (Current Slope / Open Slope)
        double slopeDenom = (Math.Abs(slopeOP) < 0.00005) ? 0.00005 : Math.Abs(slopeOP);
        double momentumRatio = Math.Abs(slopeCP) / slopeDenom;

        // C. Excess Energy (Simple Subtraction)
        double rawScore = (structureRatio - 1.0) + (momentumRatio - 1.0);

        // D. The Ghost Signal Fix (Must be positive expansion)
        if (rawScore <= 0)
            return 0.0;

        // E. Direction & Normalize
        int direction = (slopeCP >= 0) ? 1 : -1;
        double normalizedScore = Math.Tanh(rawScore * direction * 2.0);

        return normalizedScore;

        // if (normalizedScore > 0.5)
        //     return SIG.BUY;
        // else if (normalizedScore < -0.5)
        //     return SIG.SELL;

        // return SIG.HOLD;
    }

    public double efficiencyRatio(in double[] sig, int period = 14)
    {
        double net = Math.Abs(sig[SHIFT] - sig[SHIFT + period]);
        double sumAbs = 0.0;
        for (int i = SHIFT; i < SHIFT + period; i++)
            sumAbs += Math.Abs(sig[i] - sig[i + 1]);
        return (sumAbs > 0) ? net / sumAbs : 0.0;
    }

    public double vWCM_Raw(int N = 10)
    {
        double sum_force = 0.0;
        double total_vol = 0.0;
        double pipVal = _indData.PipValue;

        for (int i = SHIFT; i < N + SHIFT; i++)
        {
            double body_pips = (_indData.Close[i] - _indData.Open[i]) / pipVal;
            sum_force += body_pips * _indData.TickVolume[i];
            total_vol += _indData.TickVolume[i];
        }

        if (total_vol <= 0)
            return 0.0; // Return absolute 0 if no volume exists

        return sum_force / total_vol;
    }

    public double vWCM_Smooth(int N = 10)
    {
        double sum_force = 0.0;
        double total_vol = 0.0;
        double pipVal = _indData.PipSize;
        for (int i = SHIFT; i < N + SHIFT; i++)
        {
            double body_pips = (_indData.Close[i] - _indData.Open[i]) / pipVal;
            sum_force += body_pips * _indData.TickVolume[i];
            total_vol += _indData.TickVolume[i];
        }

        if (total_vol <= 0)
            return _indData.DBL_EPSILON;
        double raw = sum_force / total_vol;
        return Math.Tanh(raw / 10.0);
    }

    public double adxKinetic(double scale = 50.0, int shift = 1)
    {
        double adx = _indData.Adx[SHIFT];
        double normAdx = Math.Min(adx / scale, 1.0);
        return (normAdx * normAdx);
    }

    // 3. adxVector (The Compass)
    // TRUTH: "Where are we going, and do we mean it?"
    // LOGIC: Direction (+DI/-DI) weighted by Potential Energy.
    // RETURNS: -1.0 (Bearish) to +1.0 (Bullish)
    public double adxVector()
    {
        double main = _indData.Adx[SHIFT];
        double plus = _indData.AdxPlus[SHIFT];
        double minus = _indData.AdxMinus[SHIFT];

        double direction = plus - minus;
        double spread = Math.Max(Math.Abs(main - plus), Math.Abs(main - minus));

        // Use our "Smart" Potential to automatically dampen noise
        double potential = adxPotential(14);
        // Or pass 'main' if you want to optimize speed: adxPotential(main) if overloaded

        double rawSignal = direction * potential * (spread * 0.01);
        return Math.Tanh(rawSignal);
    }

    public double volatilityAnomaly()
    {
        double stdCurrent = _indData.StdClose[SHIFT];
        double avgStd = _indData.AvgStd[SHIFT];
        return (avgStd > 0) ? (stdCurrent / avgStd) : 1.0;
    }

    //+------------------------------------------------------------------+
    //| Function: Calculate MA Acceleration                              |
    //| Returns:  True if acceleration is positive (gaining strength)    |
    //+------------------------------------------------------------------+
    //bool              isTrendAccelerating(string symbol, int timeframe, int ma_period)
    public bool isTrendAccelerating(in double[] sig, int shift = 1)
    {
        // 1. Get the MA values for the last 3 completed bars
        //double ma0 = iMA(symbol, timeframe, ma_period, 0, MODE_SMA, PRICE_CLOSE, 1);
        //double ma1 = iMA(symbol, timeframe, ma_period, 0, MODE_SMA, PRICE_CLOSE, 2);
        //double ma2 = iMA(symbol, timeframe, ma_period, 0, MODE_SMA, PRICE_CLOSE, 3);

        double ma0 = sig[shift];
        double ma1 = sig[shift + 1];
        double ma2 = sig[shift + 2]; ;


        // 2. First Derivative: Velocity (Slope)
        // How much did the MA change between bars?
        double velocity_now = ma0 - ma1; // Current slope
        double velocity_prev = ma1 - ma2; // Previous slope

        // 3. Second Derivative: Acceleration
        // Is the slope getting steeper?
        double acceleration = velocity_now - velocity_prev;

        // Logic: If velocity is positive AND acceleration is positive,
        // the bullish trend is gaining strength.
        if (velocity_now > 0 && acceleration > 0)
            return (true);

        // Logic: If velocity is negative AND acceleration is negative,
        // the bearish trend is gaining strength.
        if (velocity_now < 0 && acceleration < 0)
            return (true);

        return (false);
    }


    //+------------------------------------------------------------------+
    //| FUNCTION: trendAccelStrength (Pure Physics / Injected Context)   |
    //| CHANGE: Now accepts 'atr' as a parameter.                        |
    //|         Removes internal iATR calls for speed and consistency.   |
    //+------------------------------------------------------------------+
    public double trendAccelStrength(in double[] sig, int shift = 1)
    {
        double v0 = sig[shift];
        double v1 = sig[shift + 1];
        double v2 = sig[shift + 2];

        // 1. Calculate Raw Physics (Velocity & Acceleration in Price)
        double velocity_now = v0 - v1;
        double velocity_prev = v1 - v2;
        double raw_accel = velocity_now - velocity_prev;

        // 2. Get Context (Injected Volatility)
        // Use the passed ATR. Guard against zero/negative ATR.
        double localVolatility = Math.Max(_indData.Atr[shift], _indData.Point);

        // 3. Normalize (The Universal Ratio)
        double rel_accel = raw_accel / localVolatility;

        // 4. Squash
        return Math.Tanh(rel_accel);
    }

    public double trendQuality(in double[] sig, int period = 14)
    {
        // 1. Get the Force (0-1)
        double force = adxKinetic(50.0, period);

        // 2. Get the Smoothness (0-1)
        double smooth = efficiencyRatio(sig, period);

        // 3. Combine
        // We multiply them because we need BOTH to be good.
        // If Force is high (1.0) but Smoothness is low (0.2), result is 0.2 (Bad).
        // If Smoothness is high (1.0) but Force is dead (0.1), result is 0.1 (Bad).
        return force * smooth;
    }
    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public double expansionCompressionRatio(in double fastS, in double slowS, in double SLOPEFLOOR = 0.3)
    {
        if (fastS * slowS <= 0) return 0.0;                // directional veto
        double absSlow = Math.Abs(slowS);
        if (absSlow < SLOPEFLOOR) return 0.0;              // macro-tide too flat
        if (absSlow < 0.00001) return 0.0;                 // avoid div-by-zero

        double ratio = Math.Abs(fastS) / absSlow;

        // Linear continuous mapping — Option 1 you chose
        // ratio < 0.5  → 0.00 (strong compression)
        // ratio = 1.0  → 0.50 (neutral)
        // ratio > 2.0  → 1.00 (strong expansion)
        double contrib = (ratio - 0.5) / 1.5;
        return Math.Max(0.0, Math.Min(1.0, contrib));
    }

    // Simple Physics-based 'Elasticity' calculation
    public double springForce(double currentPrice, double smaValue)
    {
        double k = 0.1; // Spring constant
        double displacement = currentPrice - smaValue;
        return -k * displacement; // Returns the 'pull' back to the mean
    }


    // Implementation
    public double slopeAccelerationRatio(
       in double fSlope,
       in double mSlope,
       in double sSlope)
    {

        double fastSlope = fSlope;
        double mediumSlope = mSlope;
        double slowSlope = sSlope;


        // Prevent division by near-zero
        if (Math.Abs(mediumSlope) < 0.000001)
            mediumSlope = 0.000001 * (mediumSlope >= 0 ? 1 : -1);
        if (Math.Abs(slowSlope) < 0.000001)
            slowSlope = 0.000001 * (slowSlope >= 0 ? 1 : -1);

        double fastOverMedium = fastSlope / mediumSlope;
        double mediumOverSlow = mediumSlope / slowSlope;

        // 1. Get the pure magnitude of the acceleration
        double ratio = Math.Abs(fastOverMedium / mediumOverSlow);

        // 2. Re-attach the directional vector (Up = Positive, Down = Negative)
        if (fastSlope < 0)
        {
            ratio = -ratio;
        }

        // Bound for stability (prevents extreme values from breaking scores)
        ratio = Math.Max(-5.0, Math.Min(5.0, ratio));

        return ratio;
    }

    //+------------------------------------------------------------------+
    //| THE GEOMETRIC FAN ENGINE: Double-Stack Acceleration              |
    //| Purpose: Validates that price is in a "Perfect Fan" expansion.   |
    //| Returns: 1.0 (Perfect Fan), 0.5 (Weak Expansion), 0.0 (Conflict) |
    //+------------------------------------------------------------------+
    public double geometricFanScore(in double fastS, in double medS, in double slowS)
    {
        // 1. THE STRUCTURAL FLOOR (Macro Tide check)
        double floor = _indData.Atr[SHIFT] * 0.30;
        double absSlow = Math.Abs(slowS);
        if (absSlow < floor) return 0.0; // The tide is too weak.

        // 2. DIRECTIONAL HARMONY (The "No Conflict" Check)
        // All three must point in the same direction.
        if (!((fastS > 0 && medS > 0 && slowS > 0) || (fastS < 0 && medS < 0 && slowS < 0)))
        {
            return 0.0;
        }

        // 3. THE RATIO STACK
        double r1 = Math.Abs(fastS) / Math.Abs(medS); // Wave vs Current
        double r2 = Math.Abs(medS) / absSlow;        // Current vs Tide

        // 4. THE VERDICT
        if (r1 >= 1.0 && r2 >= 1.0) return 1.0; // PERFECT FAN (Acceleration)
        if (r1 >= 0.8 && r2 >= 0.8) return 0.5; // COMPRESSION (Exhaustion)

        return 0.0; // GEOMETRIC DISORDER
    }

    //+------------------------------------------------------------------+
    //| FRACTAL ALIGNMENT: Grades the harmony of the moving averages     |
    //+------------------------------------------------------------------+
    public double fractalAlignment(in double fastS, in double medS, in double slowS)
    {
        // Normalize the directions to simple +1 (Up), -1 (Down), or 0 (Flat)
        // We use a microscopic threshold (0.01) to ignore completely flat MA noise
        int dFast = (fastS > 0.01) ? 1 : ((fastS < -0.01) ? -1 : 0);
        int dMed = (medS > 0.01) ? 1 : ((medS < -0.01) ? -1 : 0);
        int dSlow = (slowS > 0.01) ? 1 : ((slowS < -0.01) ? -1 : 0);

        // SCENARIO 1: Absolute Perfection (The Trend Explosion)
        // All three moving averages agree.
        if (dFast == dMed && dMed == dSlow && dFast != 0)
        {
            return 1.0;
        }

        // SCENARIO 2: The Healthy Pullback (Buy the Dip)
        // Medium and Slow agree (Macro Trend is safe), but Fast is taking a breather.
        if (dMed == dSlow && dMed != 0)
        {
            return 0.50; // We return 50% confidence. It won't veto the trade, but lot size will be smaller!
        }

        // SCENARIO 3: The Early Reversal / Macro Transition
        // Fast and Medium agree, but they are fighting the Slow Macro MA.
        // This is a riskier, early-stage breakout.
        if (dFast == dMed && dFast != 0)
        {
            return 0.25; // 25% confidence.
        }

        // SCENARIO 4: Total Chaos (Whipsaw / Choppy Market)
        return 0.0;
    }

    //+------------------------------------------------------------------+
    //| KINEMATIC ACCELERATION: Measures Trend Expansion vs Compression  |
    //+------------------------------------------------------------------+
    public double kinematicAcceleration(in double fastS, in double slowS, in double slopeFloor = 0.3)
    {
        // 1. DIRECTIONAL HARMONY
        // If they are pointing in opposite directions, it's a mess. Ratio is 0.
        if (fastS * slowS <= 0) return 0.0;

        double absSlow = Math.Abs(slowS);
        double absFast = Math.Abs(fastS);

        // 2. THE INSTITUTIONAL FLOOR
        // If the macro tide is dead, ignore it.
        if (absSlow < slopeFloor) return 0.0;

        // 3. THE ACCELERATION RATIO
        // > 1.0 means Expansion. < 1.0 means Compression.
        return (absFast / absSlow);
    }

    public double overallMarketForce(int period)
    {
        // Ensure we use SHIFT = 1 to match MQL4's 'last closed bar' logic
        int S = 1;

        double stdCurrent = _indData.StdClose[S];
        double avgStd = _indData.AvgStd[S];

        double adx = _indData.Adx[S];
        double plusDI = _indData.AdxPlus[S];
        double minusDI = _indData.AdxMinus[S];

        double rawForce = (plusDI - minusDI) * (adx / 100.0);
        double volWeight = (avgStd > 0) ? (stdCurrent / avgStd) : 1.0;

        return Math.Tanh(rawForce * volWeight * 0.1);
    }

    public double layeredMomentumFilter(in double[] values, int N = 20)
    {
        // --- Step 1: Get Smoothed Slopes ---
        double[] slopes = _stats.slopeRange_v2(values, _indData, N, 3, 1);

        int SIZE = slopes.Length;
        // if (SIZE < (N-1))
        //     return 0;

        // --- Step 2: Directional Consensus (The 80% Rule) ---
        double slopeBuy = 0;
        double slopeSell = 0;
        for (int i = 0; i < SIZE; i++)
        {
            if (slopes[i] > 0)
                slopeBuy++;
            if (slopes[i] < 0)
                slopeSell++;
        }

        SIG sig = SIG.NOTRADE;

        if (slopeBuy >= 0.8 * SIZE)
            sig = SIG.BUY;
        else if (slopeSell >= 0.8 * SIZE)
            sig = SIG.SELL;
        if (sig == SIG.NOTRADE)
            return 0;

        // --- Step 3: DYNAMIC ADX GATE (Unified Trend Power) ---
        // We reject anything below "Trend Power 0.75" (ADX 15).
        // This allows M15 scalps AND early H1 entries, but kills dead markets.
        double power = adxPotential(14);

        if (power < 0.75)
            return 0;

        // --- Step 4: Histogram Gate (Momentum Conviction) ---
        int domBin = _stats.HistogramMagnitude(slopes, N, 5, 0.2);
        if (domBin == -1)
            return 0;
        if (domBin < 2)
            return 0;

        // --- Step 5: Quality Gate (Statistical Stability) ---
        double skew = _stats.CalculateSkewness(slopes, N);
        double kurt = _stats.CalculateKurtosis(slopes, N);

        // Reject Parabolic Bubbles (High Skew > 0.5)
        if (Math.Abs(skew) > 0.5)
            return 0;
        // Reject News Spikes (High Kurtosis > 2.0)
        if (kurt > 2.0)
            return 0;

        // --- Final Signal Trigger ---
        if (sig == SIG.BUY)
            return 1.0;
        if (sig == SIG.SELL)
            return -1.0;
        return 0.0;
    }

    // 1. TIME DECAY (Linear)
    public double getLinearTimeRetention(int barsHeld, double decayRate = 0.05, double floor = 0.60)
    {
        double retention = 1.0 - (barsHeld * decayRate);
        return Math.Max(retention, floor);
    }

    // 2. VOLATILITY DECAY (Adaptive Trend Following)
    public double getVolAdaptiveRetention()
    {
        // 1. Get Normalized Volatility (0.0 to 1.0)
        //double volScore = atrStrength(atr);
        double volScore = atrKinetic();

        // 2. Trend Following Logic
        // Sqrt makes it loosen quickly as soon as volatility starts.
        double retention = 0.98 - (0.16 * Math.Sqrt(volScore));

        return Math.Max(retention, 0.70);
    }

    // 3. HYBRID DECAY (Time + Volatility)
    public double getHybridRetention(int barsHeld)
    {
        double timeRet = getLinearTimeRetention(barsHeld, 0.02, 0.80);
        double volRet = getVolAdaptiveRetention();
        return (timeRet * volRet);
    }

    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public double getHybridRetention_v2(int barsHeld, double trendQualityScore = 0.0)
    {
        // Base time decay (slower than before)
        double timeRet = 1.0 - (barsHeld * 0.015);
        timeRet = Math.Max(timeRet, 0.75);

        // Volatility adaptive (your original logic, slightly loosened)
        double volScore = atrKinetic();
        double volRet = 0.98 - (0.12 * Math.Sqrt(volScore));
        volRet = Math.Max(volRet, 0.72);

        // NEW: Trend strength bonus (Keeper mode = slower decay)
        double trendBonus = 1.0 + (trendQualityScore * 0.35);   // 0.0 → 1.35x slower decay

        return timeRet * volRet * trendBonus;
    }


    //+------------------------------------------------------------------+
    //| getPeakDecay — Unified Adaptive Peak Decay Calculator            |
    //|                                                                  |
    //| • Strategy-based: ATR, ADX, ER, or mix                           |
    //| • Output: double 0.82–0.98 (tight to loose)                     |
    //| • Use in signals: PEAK_DROP = getPeakDecay(STRAT_ATR, atr);     |
    //+------------------------------------------------------------------+
    public double getPeakDecay(DECAY_STRATEGY strat = DECAY_STRATEGY.STRAT_ATR, double period = 14, int shift = 1)
    {
        double baseDecay = 0.82;  // Min decay (tight in weak regimes)
        double scale = 0.16; // Max addition (loose in strong regimes)
        double norm = 0.0;   // Normalized strength (0–1)
        double atr = _indData.Atr[SHIFT];
        switch (strat)
        {
            case DECAY_STRATEGY.STRAT_ATR:
                {
                    // Your ATR norm (volatility fuel)

                    double pipValue = _indData.PipSize;
                    double atrPips = (pipValue > 0) ? atr / pipValue : 0.0;
                    double tfScale = (_indData._Period > 1) ? Math.Log(_indData._Period) : 1.0;
                    double atrCeiling = Math.Ceiling(12.0 * tfScale);
                    norm = Math.Min(Math.Max(atrPips / atrCeiling, 0.0), 1.0);
                    break;
                }

            case DECAY_STRATEGY.STRAT_ADX:
                {
                    // ADX norm (trend quality)
                    double adx = _indData.Adx[SHIFT];
                    norm = Math.Min(adx / 50.0, 1.0);  // 0–1 (ADX>50 rare)
                    break;
                }
            //case STRAT_ER:
            //  {
            //   // ER norm (efficiency)
            //   double net = MathAbs(close[shift] - close[shift + (int)period]);
            //   double sumAbs = 0.0;
            //   for(int i = shift; i < shift + (int)period; i++)
            //      sumAbs += MathAbs(close[i] - close[i+1]);
            //   norm = (sumAbs > 0) ? net / sumAbs : 0.0;
            //   break;
            //  }
            case DECAY_STRATEGY.STRAT_MIX:
                {
                    // Weighted mix (0.5 ATR + 0.3 ADX + 0.2 ER)
                    double atrNorm = getPeakDecay(DECAY_STRATEGY.STRAT_ATR, period, shift);
                    double adxNorm = getPeakDecay(DECAY_STRATEGY.STRAT_ADX, period, shift);
                    double erNorm = getPeakDecay(DECAY_STRATEGY.STRAT_ER, period, shift);
                    norm = 0.5 * atrNorm + 0.3 * adxNorm + 0.2 * erNorm;
                    break;
                }
        }

        // Your square-root curve (looser in strong regimes)
        double PEAK_DROP = baseDecay + scale * Math.Sqrt(norm);
        PEAK_DROP = Math.Max(Math.Min(PEAK_DROP, 0.99), 0.70);  // Clamp

        return PEAK_DROP;
    }

    //+------------------------------------------------------------------+
    //| Bayesian Next-Bar Probability (core table) – unchanged           |
    //+------------------------------------------------------------------+
    public double bayesianNextBarProb(bool volHigh, bool accelStrong)
    {
        if (volHigh && accelStrong)
            return 0.84;
        if (volHigh && !accelStrong)
            return 0.58;
        if (!volHigh && accelStrong)
            return 0.64;
        return 0.40;
    }

    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public FEATURE_VECTOR getFeatureVector(in IndData indData, in int SHIFT = 1)
    {
        FEATURE_VECTOR fV = new FEATURE_VECTOR();

        // Group A: Momentum & Velocity (The Drivers)
        //double slopeIma5 = (indData.ima5[SHIFT]-indData.ima5[SHIFT+5])/(5*pipVal);
        //double slopeIma30 = (indData.ima30[SHIFT]-indData.ima30[SHIFT+5])/(5*pipVal);
        double[] slopesArrIma5;
        double[] slopesArrIma30;
        double[] slopesArrIma60;

        slopesArrIma5 = _stats.slopeRange_v2(indData.Ima5, indData, 100, 5, SHIFT);
        slopesArrIma30 = _stats.slopeRange_v2(indData.Ima30, indData, 100, 5, SHIFT);
        slopesArrIma60 = _stats.slopeRange_v2(indData.Ima60, indData, 100, 5, SHIFT);

        //double slopeIma5 = ss.imaSlope5Data.val2;
        //double slopeIma30 = ss.imaSlope30Data.val2;


        fV.slopeIma5 = _stats.GetDistribution(slopesArrIma5, 0).zScore;
        fV.slopeIma30 = _stats.GetDistribution(slopesArrIma30, 0).zScore;
        fV.adxPlusMinusDiff = (indData.AdxPlus[SHIFT] - indData.AdxMinus[SHIFT]) / 50.0;
        fV.rsi = (indData.Rsi[SHIFT] - 50) / 50;

        // Group B: Volatility & Energy (The Fuel)
        fV.atr = _stats.GetDistribution(indData.Atr, SHIFT).zScore;
        fV.stdDevCP = indData.StdClose[SHIFT];
        fV.adx = (indData.Adx[SHIFT] / 100);
        fV.tVol = _stats.GetDistribution(indData.TickVolume, SHIFT).zScore;


        // Group C: Structure & Stretch (The Geometry)
        double[] pElastArr = new double[102];
        for (int i = 0; i < 100; i++)
        {
            pElastArr[i] = (indData.Close[i + SHIFT] - indData.Ima60[i + SHIFT]) / indData.PipSize;
        }
        //double priceElasticity = (Bid - indData.ima60[SHIFT]) / pipVal;
        fV.priceElasticity = _stats.GetDistribution(pElastArr, 0).zScore;
        //stats.zScore(pElastArr[0], stats.mean(pElastArr, 100), stats.stdDev(pElastArr, 0, 100));

        fV.mfi = (indData.Mfi[SHIFT] - 50) / 50;
        fV.vWCM = vWCM_Smooth(30);
        fV.expansionCompression = expansionCompressionRatio(slopesArrIma30[0], slopesArrIma60[0], (_indData.Atr[SHIFT] * 0.3));

        // Group D: The "Super-Feature" Additions.
        fV.bayesianScore = indData.BayesianHoldScore;
        fV.neuronScore = indData.NeuronHoldScore;
        fV.fMSR_Norm = indData.FMSR_Norm;
        fV.fractalAlignment = indData.FractalAlignment;
        return fV;
    }

    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public double marketIntensity(in FEATURE_VECTOR fV)
    {

        double[] cloud16D = new double[16]; // Declare the container

        // Group A: Momentum
        cloud16D[0] = fV.slopeIma5;
        cloud16D[1] = fV.slopeIma30;
        cloud16D[2] = fV.adxPlusMinusDiff;
        cloud16D[3] = fV.rsi;

        // Group B: Energy
        cloud16D[4] = fV.atr;
        cloud16D[5] = fV.stdDevCP;
        cloud16D[6] = fV.adx;
        cloud16D[7] = fV.tVol;

        // Group C: Geometry
        cloud16D[8] = fV.priceElasticity;
        cloud16D[9] = fV.mfi;
        cloud16D[10] = fV.vWCM;
        cloud16D[11] = fV.expansionCompression;

        // Group D: Super-Features
        cloud16D[12] = fV.bayesianScore;
        cloud16D[13] = fV.neuronScore;
        cloud16D[14] = fV.fMSR_Norm;
        cloud16D[15] = fV.fractalAlignment;

        double globalIntensity = Vector<double>.Build.Dense(cloud16D).L2Norm();
        //double globalIntensity = cloud16D.L2Norm();
        return globalIntensity;
    }
    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public double marketRegime(in FEATURE_VECTOR fV)
    {

        // --- 1. COORDINATE CALCULATION (The Position) ---
        double posX = (fV.slopeIma5 + fV.slopeIma30 + fV.adxPlusMinusDiff + fV.rsi + fV.fractalAlignment);
        double posY = (fV.atr + fV.stdDevCP + fV.adx + fV.tVol);
        double posZ = (fV.priceElasticity + fV.mfi + fV.vWCM + fV.expansionCompression);

        // 2. The 3D Projected Norm
        double[] projection3D = new double[3];
        projection3D[0] = posX; // These were already calculated as sums
        projection3D[1] = posY;
        projection3D[2] = posZ;
        // System.Console.WriteLine("[ [MOMENTUM] - X:"+ posX:F2| [ENERGY] - Y: {posY:F2} | [STRUCTURE] - Z: {posZ:F2} ]");
        double regimeMagnitude = Vector<double>.Build.Dense(projection3D).L2Norm();
        return regimeMagnitude;
    }


    //+------------------------------------------------------------------+
    //|                                                                  |
    //+------------------------------------------------------------------+
    public double bayesianHoldScore(
   in double[] maArray, in double[] closeArray, in double[] openArray, in double[] volumeArray,
   int barsHeld, double atr, bool useOverallForce = true)
    {
        // 1. Gather Atoms
        double tq = trendQuality(maArray);

        // OLD (Bullish only)
        // bool accelStrong = (trendAccelStrength(maArray) > 0.42);

        // 2. Fix Acceleration Directionality
        // OLD: bool accelStrong = (trendAccelStrength(...) > 0.42);
        // NEW: Check Magnitude.
        bool accelStrong = (Math.Abs(trendAccelStrength(maArray, SHIFT)) > 0.10);

        double volRaw = vWCM_Raw(10);
        // 1. Fix Volume Directionality
        // OLD: bool volHigh = (volRaw > 0.35);  <-- Fails on Sells
        // NEW: Check Magnitude. Also lowered threshold slightly for M15.
        bool volHigh = (Math.Abs(volRaw) > 0.15);


        double anomaly = volatilityAnomaly();
        // Soft anomaly penalty: Only punish if > 1.5 (extreme shock)
        double anomalyMultiplier = (anomaly > 1.5) ? 0.7 : 1.0;

        double regimeMultiplier = 1.0;
        if (useOverallForce)
        {
            double omf = overallMarketForce(14);
            regimeMultiplier = (omf > 0.35) ? 1.2 : (omf < -0.35 ? 0.8 : 1.0);
        }

        // 2. Likelihood Ratios (Evidence Weight)
        // LR = 1.0 (Neutral)
        double lr = 1.0;
        // Trend Quality: The strongest predictor
        // Adjusted Likelihood Ratios in bayesianHoldScore
        // Trend Quality Logic
        if (tq > 0.20)
            lr *= 2.0; // Unicorn
        else if (tq > 0.12)
            lr *= 1.2; // <-- ADD THIS: The "Grinder" Boost. (0.16 falls here now)
        else if (tq < 0.08)
            lr *= 0.4; // Swamp
        else
            lr *= 0.9; // Weak Drift

        // Acceleration: Good for short term
        lr *= (accelStrong ? 1.4 : 0.8);

        // Volume: Confirmation
        lr *= (volHigh ? 1.2 : 0.9);

        // Penalties
        lr *= anomalyMultiplier;
        lr *= regimeMultiplier;

        // 3. Time Decay (The Prior eroding over time)
        // We reduce the LR effective power by 2% per bar
        double timeDecay = Math.Max(0.2, 1.0 - (barsHeld * 0.02));
        lr *= timeDecay;

        // 4. Probability Conversion (Sigmoid-like via Odds)
        // P = Odds / (1 + Odds)

        //Print("Bayesian Hold | AccelRaw=", DoubleToStr(trendAccelStrength(maArray,atr), 4),
        //      " | TQRaw=", DoubleToStr(tq, 4),
        //      " | VolRaw=", DoubleToStr(volRaw, 4));

        return lr / (1.0 + lr);
    }




    //+------------------------------------------------------------------+
    //| NEURON HOLD SCORE — ALIGNED WITH CONTINUOUS fMSR                |
    //| FIX: Uses closed bars [1], continuous fMSR, and rewards volume  |
    //| intensity correctly.                                            |
    //+------------------------------------------------------------------+
    public double neuronHoldScore(
       in double[] maArray, in double[] closeArray, in double[] openArray, in double[] volumeArray,
       int barsHeld, double atr, bool useOverallForce = true)
    {
        double tq = trendQuality(maArray);

        // 1. Universal Acceleration (Magnitude)
        bool accelStrong = (Math.Abs(trendAccelStrength(maArray)) > 0.18);

        // 2. Universal Volume (Magnitude) — your good fix
        double volRaw = vWCM_Raw(10);
        bool volHigh = (Math.Abs(volRaw) > 0.15);

        double bayesP = bayesianNextBarProb(volHigh, accelStrong);
        double retention = getHybridRetention_v2(barsHeld, tq);
        double anomaly = volatilityAnomaly();
        double anomalyFactor = (anomaly > 1.0) ? Math.Max(0.82, 1.0 - (anomaly - 1.0) * 0.22) : 1.0;

        double regimeBias = 1.0;
        if (useOverallForce)
        {
            double omf = overallMarketForce(14);
            regimeBias = (omf > 0.35) ? 1.15 : (omf < -0.35 ? 0.75 : 1.0);
        }

        // === CONTINUOUS fMSR (Aligned with RefreshPhysicsData) ===
        // Use closed bar [1] to avoid repainting on M15 USDJPY
        double pipVal = _indData.PipSize;
        double fastSlope = (maArray[1] - maArray[4]) / (3 * pipVal);
        double medSlope = (maArray[1] - maArray[11]) / (10 * pipVal);
        double slowSlope = (maArray[1] - maArray[31]) / (30 * pipVal);  // consistent with RefreshPhysicsData

        double fMSR_raw = slopeAccelerationRatio(fastSlope, medSlope, slowSlope);
        double abs_fMSR = Math.Abs(fMSR_raw); // <-- ADD THIS
                                              // Continuous linear version — Option 1 (smooth 0.0 → 1.0)
                                              //   double fMSR_norm = MathMax(0.0, MathMin(1.0, (fMSR_raw - 0.5) / 1.5));
        double fMSR_norm = 0.0;
        // Bimodal Activation (U-Shape)
        if (abs_fMSR >= 0.40)
        {
            // Perfect Expansion
            fMSR_norm = 1.0;
        }
        else if (abs_fMSR <= 0.15)
        {
            // Perfect Compression (Energy Storage)
            // We give the Neural Network full points because a squeeze is a highly valid setup!
            fMSR_norm = 1.0;
        }
        else
        {
            // The Toxic No-Go Zone
            // Penalize the Neural score so it doesn't fire in choppy markets
            fMSR_norm = 0.0;
        }

        // ===============================================
        // SCORE COMBINATION
        // ===============================================
        double score = 0.0;

        score += tq * 0.40;
        score += bayesP * 0.25;
        score += (retention - 0.7) / 0.7 * 0.12;
        score += anomalyFactor * 0.10;
        score += (regimeBias - 0.9) * 0.08;
        score += fMSR_norm * 0.15;          // ← Now active and continuous
        score += Math.Abs(volRaw) * 0.10;    // Your Sell Volume Intensity fix

        // Sigmoid squash (smooth 0.0–1.0 output)
        score = 1.0 / (1.0 + Math.Exp(-10.0 * (score - 0.5)));

        return Math.Max(Math.Min(score, 1.0), 0.0);
    }


    //+------------------------------------------------------------------+
    //| THE HYPERBOLIC MAP (Bimodal Squeeze Authorization)               |
    //+------------------------------------------------------------------+
    public int getHyperbolicCombinedScore(double b, double n, double f, double fra)
    {

        const double TRADE_EXPANSION = 0.40;
        const double TRADE_COMPRESSION = 0.15;

        double absF = Math.Abs(f);

        // 1. THE GUARD CLAUSE (No-Go Zone)
        if ((absF > TRADE_COMPRESSION) && (absF < TRADE_EXPANSION)) return 0;

        // --- PHASE 1: EXPANSION (The Trend) ---
        if (absF >= TRADE_EXPANSION)
        {
            // Additive structure: We add probabilities AND fractal alignment
            double combined = (b * 1.5) + (n * 1.2) + (fra * 1.0);

            // We subtract an offset (e.g., 2.0) so the tanh only goes positive
            // if the combined score is extremely high.
            double score = Math.Tanh(combined - 2.0);

            if (score > 0.20) return 1;   // Authorize
            if (score < -0.20) return -1; // Veto
            return 0;
        }

        // --- PHASE 2: COMPRESSION (The Squeeze) ---
        else if (absF <= TRADE_COMPRESSION)
        {
            // Ignore Fractals. We boost the weight of the Smart Money probabilities (b & n)
            double combined = (b * 1.8) + (n * 1.5);

            double score = Math.Tanh(combined - 1.8);

            // We require a slightly higher tanh output to authorize a squeeze in the dark
            if (score > 0.30) return 1;  // Authorize Vanguard Snipe
            if (score < -0.30) return -1; // Veto toxic compression
            return 0;
        }

        return 0;
    }

    //+------------------------------------------------------------------+
    //| THE COBB-DOUGLAS MAP (Dual-Phase Squeeze Authorization)          |
    //+------------------------------------------------------------------+
    public int getCobbDouglasCombinedScore(double b, double n, double f, double fra)
    {

        // 1. MUST BE DOUBLE.
        const double TRADE_EXPANSION = 0.4;
        // Allow for standard M15 MA noise to be classified as a squeeze
        const double TRADE_COMPRESSION = 0.15;

        // 2. SAFEGUARD NEGATIVE TRENDS
        double absF = Math.Abs(f);

        // 3. THE GUARD CLAUSE (No-Go Zone)
        if ((absF > TRADE_COMPRESSION) && (absF < TRADE_EXPANSION)) return 0;

        // --- PHASE 1: EXPANSION (The Trend) ---
        if (absF >= TRADE_EXPANSION)
        {
            double trendConf = Math.Pow(n + 0.01, 1.2) * Math.Pow(b + 0.01, 1.5) * (fra + 0.01);

            if (trendConf > 0.20) return 1;
            if (trendConf < 0.05) return -1;
            return 0;
        }

        // --- PHASE 2: COMPRESSION (The Squeeze) ---
        else if (absF <= TRADE_COMPRESSION)
        {
            double squeezeConf = Math.Pow(n + 0.01, 1.2) * Math.Pow(b + 0.01, 1.5);

            if (squeezeConf > 0.35)
            {
                return 1;
            }
            if (squeezeConf < 0.15) return -1;
            return 0;
        }

        return 0;
    }

}
