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


public enum DECAY_STRATEGY
{
    STRAT_ATR,       // Volatility fuel (your current)
    STRAT_ADX,       // Trend quality
    STRAT_ER,        // Efficiency
    STRAT_MIX        // Weighted mix (0.5 ATR + 0.3 ADX + 0.2 ER)
}




public struct DTYPE
{
    public double val1;
    public double val2;
    public double val3;
    public double val4;
    public double val5;
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

}

public struct STRATEGY_STATE
{
    public bool newCandle { get; set; }
    public bool inTrade { get; set; }
}