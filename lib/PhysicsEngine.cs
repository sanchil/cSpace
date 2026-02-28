using System;
using System.Net.NetworkInformation;

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
    public double adxPotential(int period = 14);

    public double VolatilityEfficiency();

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
    public double FMSR { get; init; }
    public double FractalAlignment { get; init; }
    public double MicroLots { get; init; }
    public double ConvictionFactor { get; init; }
    public int PhysicsAction { get; init; }
    public double PipValue { get; init; }
    public double Current_Period { get; init; }
}

public record FeatureVector()
{
    public double Rsi { get; init; }
    public double Atr { get; init; }
    public double BayesianScore { get; init; }
}

public class PhysicsEngine : IPhysicsEngine
{
    private IndData _indData;
    private Stats _stats;
    private Utils _utils;

    private int SHIFT = 0; // Default shift for indicator calculations

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
        double pipValue = _indData.PipValue;
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

}
