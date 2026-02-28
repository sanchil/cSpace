using System.Reflection.Metadata;

namespace Phy.Lib;

public interface IStats
{
    public double[] slopesVal(
        in double[] sig,
        in int SLOPEDENOM = 3,
        in int SLOPEDENOM_WIDE = 5,
        in double pipValue = 0.0001,
        in int shift = 1);

}

public class Stats : IStats
{
    public double[] slopesVal(
            in double[] sig,
            in int SLOPEDENOM = 3,
            in int SLOPEDENOM_WIDE = 5,
            in double pipValue = 0.0001,
            in int shift = 0)
    {
        double[] slopes = new double[4];
        double ppip = pipValue;
        if (ppip == 0.0) ppip = 0.0001; // Default to a standard pip value if not provided
        slopes[0] = (sig[shift] - sig[shift + SLOPEDENOM]) / (SLOPEDENOM * ppip);
        slopes[1] = (sig[shift] - sig[shift + SLOPEDENOM_WIDE]) / (SLOPEDENOM_WIDE * ppip);
        slopes[2] = (sig[shift] - sig[sig.Length - 1]) / ((sig.Length - 1) * ppip);
        slopes[3] = slopes[0] - slopes[1]; // Slope difference as a measure of momentum change


        return slopes;
    }


}