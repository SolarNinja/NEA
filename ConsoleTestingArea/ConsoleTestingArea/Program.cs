using System.Numerics;

namespace ConsoleTestingArea
{
    internal class Program
    {
        static void Main(string[] args)
        {
            Tree tree = new Tree(new FixedBody("Sol", new Vector3(Convert.ToSingle(2.5 * Math.Pow(10, 7)), Convert.ToSingle(2.5 * Math.Pow(10, 7)), 0), 5.972 * Math.Pow(10, 24)));

            List<string> heritage = new List<string>();
            heritage.Add("Sol");

            Vector3 calculatedAbsolutePosition = new Vector3();
            bool failed = false;
            tree.AddToTree(heritage, new MovingBody(ref failed, tree, "Earth", 1, 1, 1, 1, new Vector3(1000, 1000, 1000), new Vector3(1, 1, 1), calculatedAbsolutePosition));

            tree.AddToTree(heritage, new MovingBody(ref failed, tree, "Mars", 1, 1, 1, 1, new Vector3(2000, 2000, 2000), new Vector3(1, 1, 1), calculatedAbsolutePosition));

            heritage.Add("Earth");

            tree.AddToTree(heritage, new MovingBody(ref failed, tree, "Moon", 1, 1, 1, 1, new Vector3(1000, -1000, 1000), new Vector3(1, 1, 1), calculatedAbsolutePosition));

            heritage.Add("Moon");

            tree.AddToTree(heritage, new MovingBody(ref failed, tree, "Moon Moon", 1, 1, 1, 1, new Vector3(10000, -10000, 10000), new Vector3(1, 1, 1), calculatedAbsolutePosition));

            heritage.Remove("Moon");
            heritage.Remove("Earth");
            heritage.Add("Mars");

            tree.AddToTree(heritage, new MovingBody(ref failed, tree, "Moon Moon Moon", 1, 1, 1, 1, new Vector3(10000, -10000, 10000), new Vector3(1, 1, 1), calculatedAbsolutePosition));

            tree.AreObjectsWithinRadius(3, new Vector3(2, 2, 2));
        }
    }
}