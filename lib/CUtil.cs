namespace Phy.Lib;

public interface IUtils
{

    public bool IsNewBar(DateTime currentBarTime);
}
public class CUtils : IUtils
{
    private DateTime _lastBar;

    public CUtils()
    {
        _lastBar = DateTime.MinValue;
    }

    public bool IsNewBar(DateTime currentBarTime)
{
    if (currentBarTime != _lastBar)
    {
        _lastBar = currentBarTime;
        return true;
    }
    return false;
}

}