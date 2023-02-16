using System.Xml;

namespace Test
{
    internal class Program
    {
        static void Main(string[] args)
        {
            /*double midpointDistance = 149599969247.9592;
            double V = -443403500.23771751;
            double A = -4454316460780803;
            double timeOffset = 0;
            double startPositionTimeOffset = 0;
            double K = 5.2669839591951317E+36;
            double G = 6.67 * Math.Pow(10, -11);
            double orbitalPeriod = 31564347.516208;
            double ParentMass = 1.9890000000000002E+30;
            double midpointTimeAdjustment = 15782173.758104002;


            string[] outputs = new string[10001];
            for (double i = 1.47009 * Math.Pow(10, 11); i < 1.471009 * Math.Pow(10, 11); i += 2 * 1000)
            {
                double timePositive = double.NaN;
                double timeNegative = double.NaN;

                double U = 2 * V * Math.Pow(i, 2) + 2 * G * ParentMass * i - Math.Pow(A, 2);

                if (i < midpointDistance)
                {
                    timePositive = (-G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + timeOffset + orbitalPeriod - startPositionTimeOffset) % orbitalPeriod;
                }
                else if (i > midpointDistance)
                {
                    timePositive = (G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + midpointTimeAdjustment + timeOffset + orbitalPeriod - startPositionTimeOffset) % orbitalPeriod;
                }
                
                if (i < midpointDistance)
                {
                    timeNegative = orbitalPeriod - (-G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + timeOffset + orbitalPeriod + startPositionTimeOffset) % orbitalPeriod;
                }
                else if (i > midpointDistance)
                {
                    timeNegative = orbitalPeriod - (G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + midpointTimeAdjustment + timeOffset + orbitalPeriod + startPositionTimeOffset) % orbitalPeriod;
                }

                outputs[Convert.ToInt64(i - 1.47009 * Math.Pow(10, 11))/10000] = Convert.ToString(timePositive);
                outputs[Convert.ToInt64(i- 1.47009 * Math.Pow(10, 11))/10000 + 1] = Convert.ToString(timeNegative);
                //outputs[Convert.ToInt64(i - 1.47009 * Math.Pow(10, 11)) / 10000] = Convert.ToString(i-1.47 * Math.Pow(10, 11));
                //outputs[Convert.ToInt64(i - 1.47009 * Math.Pow(10, 11)) / 10000 + 1] = Convert.ToString(i-1.47 * Math.Pow(10, 11));
            }

            File.WriteAllLines("output.txt", outputs);*/

            Console.WriteLine(ulong.MaxValue / (60 * 60 * 24 * 365));
        }
    }
}