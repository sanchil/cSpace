using System;
using System.Runtime.Serialization;
using cAlgo.API;
using cAlgo.API.Indicators;
using Phy.Lib;

namespace Phy.Bot
{
    [Robot(AccessRights = AccessRights.None)]
    public class PhyBot : Robot, IBotEngine
    {
        private IndData _indData;
        private PhysicsEngine _engine;
        private Strategy _strategy;

        private Stats _stats;
        private Utils _utils;

        // Results objects for indicators
        private StandardDeviation _stdClose;
        private StandardDeviation _stdOpen;
        private MoneyFlowIndex _mfi;
        private AverageTrueRange _atr;
        private AverageDirectionalMovementIndexRating _adx;
        private OnBalanceVolume _obv;
        private RelativeStrengthIndex _rsi;
        private SimpleMovingAverage _ima5;
        private SimpleMovingAverage _ima14;
        private SimpleMovingAverage _ima30;
        private SimpleMovingAverage _ima60;
        private SimpleMovingAverage _ima120;
        private SimpleMovingAverage _ima240;
        private SimpleMovingAverage _ima500;


        public double[] GetHistory(Func<int, double> getPricePtr, int count, int shift = 1)
        {
            double[] buffer = new double[count];
            int startPoint = Bars.Count - 1 - shift;

            for (int i = 0; i < count; i++)
            {
                int targetIndex = startPoint - i;
                if (targetIndex >= 0)
                {
                    double val = getPricePtr(targetIndex);

                    // If the indicator is still "warming up", default to 0.0
                    buffer[i] = double.IsNaN(val) ? 0.0 : val;
                }
                else
                {
                    buffer[i] = 0.0;
                }
            }
            return buffer;
        }
        // FIX: Added 'void' return type to match the interface contract
        public void InitIndData()
        {
            _indData = new IndData
            {
                MagicNumber = 1002,
                BarsHeld = 0,
                Open = GetHistory(idx => Bars.OpenPrices[idx], 500),
                High = GetHistory(idx => Bars.HighPrices[idx], 120),
                Low = GetHistory(idx => Bars.LowPrices[idx], 120),
                Close = GetHistory(idx => Bars.ClosePrices[idx], 500),
                Time = GetHistory(idx => Bars.OpenTimes[idx].ToOADate(), 500).Select(t => DateTime.FromOADate(t)).ToArray(),
                TickVolume = GetHistory(idx => Bars.TickVolumes[idx], 500),
                StdClose = GetHistory(idx => _stdClose.Result[idx], 500),
                StdOpen = GetHistory(idx => _stdOpen.Result[idx], 500),
                Mfi = GetHistory(idx => _mfi.Result[idx], 500),
                Atr = GetHistory(idx => _atr.Result[idx], 500),
                Adx = GetHistory(idx => _adx.ADXR[idx], 500),     // main line = ADXR (smoothed ADX)
                AdxPlus = GetHistory(idx => _adx.DIPlus[idx], 500),   // +DI
                AdxMinus = GetHistory(idx => _adx.DIMinus[idx], 500),
                Obv = GetHistory(idx => _obv.Result[idx], 500),
                Rsi = GetHistory(idx => _rsi.Result[idx], 500), // -DI
                Ima5 = GetHistory(idx => _ima5.Result[idx], 500),
                Ima14 = GetHistory(idx => _ima14.Result[idx], 500),
                Ima30 = GetHistory(idx => _ima30.Result[idx], 500),
                Ima60 = GetHistory(idx => _ima60.Result[idx], 500),
                Ima120 = GetHistory(idx => _ima120.Result[idx], 500),
                Ima240 = GetHistory(idx => _ima240.Result[idx], 500),
                Ima500 = GetHistory(idx => _ima500.Result[idx], 500),
                PipValue = Symbol.PipValue,
                Current_Period = 60.0,
                Shift = 0
            };
        }

        public void printData(in IndData data)
        {
            Print("=== IndData Snapshot ===");
            Print($"MagicNumber: {data.MagicNumber}");
            Print($"Open: {data.Open[0]}");
            Print($"Close: {data.Close[0]}");
            Print($"StdClose: {data.StdClose[0]}");
            Print($"StdOpen: {data.StdOpen[0]}");
            Print($"Mfi: {data.Mfi[0]}");
            Print($"Atr: {data.Atr[0]}");
            Print($"Adx: {data.Adx[0]}");
            Print($"AdxPlus: {data.AdxPlus[0]}");
            Print($"AdxMinus: {data.AdxMinus[0]}");
            Print($"Obv: {data.Obv[0]}");
            Print($"Rsi: {data.Rsi[0]}");
            Print($"Ima5: {data.Ima5[0]}");
            Print($"Ima14: {data.Ima14[0]}");
            Print($"Ima30: {data.Ima30[0]}");
            Print($"Ima60: {data.Ima60[0]}");
            Print($"Ima120: {data.Ima120[0]}");
            Print($"Ima240: {data.Ima240[0]}");
            Print($"Ima500: {data.Ima500[0]}");
            Print($"TradePosition: {data.TradePosition}");
            Print($"CurrSpread: {data.CurrSpread}");
            Print($"Shift: {data.Shift}");
            Print($"BarsHeld: {data.BarsHeld}");
            Print($"HoldScore: {data.HoldScore}");

            // Optionally print the latest values of the historical arrays
            if (data.Close.Length > 0)
                Print($"Latest Close Price: {data.Close[0]}");
            if (data.StdClose.Length > 0)
                Print($"Latest StdDev of Close: {data.StdClose[0]}");
            if (data.Mfi.Length > 0)
                Print($"Latest MFI: {data.Mfi[0]}");
        }

        protected override void OnStart()
        {
            Print("Phy.Bot initialized. Connecting to Physics Engine...");
            _stdClose = Indicators.StandardDeviation(Bars.ClosePrices, 20, MovingAverageType.Simple);
            _stdOpen = Indicators.StandardDeviation(Bars.OpenPrices, 20, MovingAverageType.Simple);
            _mfi = Indicators.MoneyFlowIndex(20);
            _atr = Indicators.AverageTrueRange(20, MovingAverageType.Simple);
            _adx = Indicators.AverageDirectionalMovementIndexRating(20);
            _obv = Indicators.OnBalanceVolume(Bars.ClosePrices);
            _rsi = Indicators.RelativeStrengthIndex(Bars.ClosePrices, 14);
            _ima5 = Indicators.SimpleMovingAverage(Bars.ClosePrices, 5);
            _ima14 = Indicators.SimpleMovingAverage(Bars.ClosePrices, 14);
            _ima30 = Indicators.SimpleMovingAverage(Bars.ClosePrices, 30);
            _ima60 = Indicators.SimpleMovingAverage(Bars.ClosePrices, 60);
            _ima120 = Indicators.SimpleMovingAverage(Bars.ClosePrices, 120);
            _ima240 = Indicators.SimpleMovingAverage(Bars.ClosePrices, 240);
            _ima500 = Indicators.SimpleMovingAverage(Bars.ClosePrices, 500);

            InitIndData();

            _utils = new Utils();
            _stats = new Stats();
            _engine = new PhysicsEngine(_indData, _stats, _utils);
            _strategy = new Strategy(_engine, _stats, _utils);

        }

        protected override void OnTick()
        {
            // High-frequency math goes here

        }

        protected override void OnBar()
        {
            InitIndData();
            _engine.SetIndData(_indData);
            onBarTask1();

        }

        protected override void OnStop()
        {
            Print("Cleaning up resources...");
        }

        private int GetBarAge(Position position)
        {
            // 1. Get the index of the bar when the position was opened
            int entryBarIndex = Bars.OpenTimes.GetIndexByTime(position.EntryTime);

            // 2. The current bar index is always the last one in the series
            int currentBarIndex = Bars.Count - 1;

            // 3. The difference is the number of bars the trade has existed
            return currentBarIndex - entryBarIndex;
        }
        void onBarTask1()
        {
            SIG signal = _strategy.GetSignal();

            Print($"SIG: {signal}");

            // Define your volume (Example: 10,000 units = 0.10 lots)
            double volumeUnits = Symbol.QuantityToVolumeInUnits(0.1);
            string label = "PhyBot_Signal";

            int allPositions = Positions.Count;

            // Count ONLY positions opened by this bot (using your "PhyLabel")
            var botPositions = Positions.FindAll(label, SymbolName);
            int activeTradesCount = botPositions.Length;
            //################## CLOSE LOGIC ##################
            if (activeTradesCount > 0)
            {
                foreach (var pos in botPositions)
                {
                    _indData = _indData with { TradePosition = (pos.TradeType == TradeType.Buy) ? SIG.BUY : SIG.SELL, BarsHeld = GetBarAge(pos) };
                    _engine.SetIndData(_indData);
                   // printData(_indData);

                    // Check if we need to close the position based on the new signal
                    if (pos.TradeType == TradeType.Buy && signal == SIG.SELL)
                    {
                        Print(">>> Closing BUY position due to SELL signal...");
                        ClosePosition(pos);
                    }
                    else if (pos.TradeType == TradeType.Sell && signal == SIG.BUY)
                    {
                        Print(">>> Closing SELL position due to BUY signal...");
                        ClosePosition(pos);
                    }
                }

            }
            //###########################################################

            //################## OPEN LOGIC #############################

            // if (activeTradesCount == 0)
            // {
            switch (signal)
            {
                case SIG.BUY:
                    Print(">>> BUY signal generated!");
                    // ExecuteMarketOrder(TradeType, SymbolName, Volume, Label, StopLoss, TakeProfit)
                    // ExecuteMarketOrder(TradeType.Buy, SymbolName, volumeUnits, label, 10, 20);
                    ExecuteMarketOrder(TradeType.Buy, SymbolName, volumeUnits, label, null, null);
                    break;

                case SIG.SELL:
                    Print(">>> SELL signal generated!");
                    // ExecuteMarketOrder(TradeType.Sell, SymbolName, volumeUnits, label, 10, 20);
                    ExecuteMarketOrder(TradeType.Sell, SymbolName, volumeUnits, label, null, null);
                    break;

                case SIG.HOLD:
                    // Do nothing, or monitor existing positions
                    break;
            }
            //           }

        }
    }
}