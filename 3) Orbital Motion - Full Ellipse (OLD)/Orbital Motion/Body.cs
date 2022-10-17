using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Windows;
using System.Reflection.Metadata.Ecma335;
using System.DirectoryServices;
using System.Linq;
using System.Windows.Forms;

namespace Orbital_Motion
{
    internal class Body
    {
        private Vector2 position;
        private Vector2 velocity;
        private Vector2 radialPosition;

        private double startingDistance;

        // Constants
        double Mass;
        static double TotalEnergy, AngularMomentum;
        //static double G = 6.6743 * Math.Pow(10, -11);
        static double G = 6.67 * Math.Pow(10, -11);

        // Psuedo constants (constant unless a change in orbit occurs)
        private double eccentricity, semiMajorAxis, semiMinorAxis, periapsis, apoapsis, orbitalPeriod;
        private Dictionary<double, double[]> predeterminedPoints;

        public Body(Vector2 position, Vector2 velocity, double Mass, double[] parentBody)
        {
            double kineticEnergy, gravitationalPotentialEnergy;

            this.position = position;
            this.velocity = velocity;
            this.Mass = Mass;

            startingDistance = Math.Sqrt(Math.Pow(position.X - parentBody[0], 2) + Math.Pow(position.Y - parentBody[1], 2));

            gravitationalPotentialEnergy = -1 * G * parentBody[2] * Mass / startingDistance;
            kineticEnergy = 0.5 * Mass * Vector2.Dot(velocity, velocity);
            // Ek = -1/2 Ep, and hence Et = 1/2 Ep, but unfortunately computer no like so I gotta do differently.

            TotalEnergy = gravitationalPotentialEnergy + kineticEnergy;
            //TotalEnergy = 1/2 * gravitationalPotentialEnergy;
            AngularMomentum = Mass * velocity.Length() * startingDistance;

            eccentricity = Math.Sqrt(1 + 2 * TotalEnergy * Math.Pow(AngularMomentum, 2) / (Math.Pow(Mass, 3) * Math.Pow(G * parentBody[2], 2)));
            semiMajorAxis = 1 / (2 / startingDistance - Vector2.Dot(velocity,velocity) / (G * parentBody[2]));
            semiMinorAxis = semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity, 2));
            // These are minorly adjusted to make sure no errors occur when calculating the extreme cases.
            periapsis = semiMajorAxis * (1 - eccentricity) + 0.01;
            apoapsis = semiMajorAxis * (1 + eccentricity) - 0.01;

            // Orbital period is defined here due to being more convienient in timing with the maths.
            predeterminedPoints = createPoints(10000, parentBody, periapsis, apoapsis, 0, out orbitalPeriod);
        }

        public void tickMovement(double secondsPerTick, Vector2 acceleration)
        {
            Vector2 acceleratedVelocity = acceleration;

            // Update approximation as currently incorrect.
            position += Vector2.Multiply(velocity, Convert.ToSingle(secondsPerTick)) + Vector2.Multiply(acceleratedVelocity, Convert.ToSingle(0.5 * secondsPerTick * secondsPerTick));
            velocity += acceleratedVelocity;
        }

        public double createConic(double[] parentBody)
        {
            // The standard gravitation parameter, -GM.
            double mew = -1 * G * parentBody[2];

            // The angular momentum.
            double l = Mass * (velocity.X * position.Y - velocity.Y * position.X);

            // The eccentricity.
            double e = Math.Sqrt(1 + (2 * TotalEnergy * l * l) / (Mass * Mass * Mass * mew * mew));

            // Interestingly, the eccentricity vector of a craft of position p and velocity v is given by a different equation.
            Vector2 eVector = Vector2.Divide(velocity * velocity * position, Convert.ToSingle(mew)) - Vector2.Divide(position, Convert.ToSingle(position.Length()));

            // The true anomaly, as calculated using the eccentricity vector.
            double theta = Math.Acos(Vector2.Dot(eVector, position) / (eVector.Length() * position.Length()));

            // The range between the body and the parent.
            double r = (l * l) / (Mass * Mass * mew) * 1 / (1 + e * Math.Cos(theta));

            return r;
        }

        // This method weights distance equally not time, so may be a little janky further away from the planet given more time spent at the "same" place.
        private Dictionary<double, double[]> createPoints(int numberOfIntervals, double[] parentBody, double periapsis, double apoapsis, double startTime, out double orbitalPeriod)
        {
            Dictionary<double, double[]> predeterminedPoints = new Dictionary<double, double[]>();

            double V = TotalEnergy / Mass;
            double A = AngularMomentum / Mass;
            double K = Math.Pow(G * parentBody[2], 2) + 2 * V * Math.Pow(A, 2);

            double midpointR = -G * parentBody[2] / (V * 2) + 0.01;
            double midpointAdjustedU = 2 * V * Math.Pow(midpointR, 2) + 2 * G * parentBody[2] * midpointR - Math.Pow(A, 2);
            double periapsisU = 2 * V * Math.Pow(periapsis, 2) + 2 * G * parentBody[2] * periapsis - Math.Pow(A, 2);
            double partA = G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * midpointAdjustedU / K), 5)));
            double midpointTimeAdjustment = -2 * partA;
            orbitalPeriod = 2 * (G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * periapsisU / K), 5))) + 1 / (2 * V) * Math.Sqrt(periapsisU) + midpointTimeAdjustment);

            double range = apoapsis - periapsis;

            // Doesnt work at exact half.
            for (int i = 0; i < 2 * numberOfIntervals; i += 2)
            {
                double R = periapsis + i * range / (2 * numberOfIntervals);
                double U = 2 * V * Math.Pow(R, 2) + 2 * G * parentBody[2] * R - Math.Pow(A, 2);
                double time = 0;

                // I actually believe this to be more efficient than defining repeating parts as it saves on read-writes.
                // Also note the use of the round; this is because an error would occur at halfway where, due to precision errors, the arcsin would recieve a value ever so slightly above 1. I might need to revisit this error to make it's handling more efficient.
                if (R < midpointR)
                {
                    time = -1 * G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K),5))) + 1 / (2 * V) * Math.Sqrt(U);
                }
                else if (R >= midpointR)
                {
                    time = G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + midpointTimeAdjustment;
                }

                predeterminedPoints.Add(i, new double[] {R, time});
                predeterminedPoints.Add(i + 1, new double[] { R, orbitalPeriod - time });
            }

            return predeterminedPoints;
        }

        public void returnValues(double time, double[] parentBody, out Vector2 radialPositionOut, out Vector2 relativeVelocity, out double angularVelocity)
        {
            double predeterminedMethod()
            {
                double[] predeterminedAdjustedTimes = new double[predeterminedPoints.Count];

                time %= orbitalPeriod;

                int count = 0;
                foreach (var point in predeterminedPoints)
                {
                    predeterminedAdjustedTimes[count] = Math.Abs(predeterminedPoints[count][1] - time);
                    count++;
                }

                double[] predeterminedPoint = predeterminedPoints[Array.IndexOf(predeterminedAdjustedTimes, predeterminedAdjustedTimes.Min())];
                return predeterminedPoint[0];
            }

            double distance = predeterminedMethod();

            double angle;

            // Like before, I actually believe this to be more efficient than defining repeating parts as it saves on read-writes.
            if (time < orbitalPeriod / 2)
            {
                angle = Math.Acos((semiMajorAxis * (1 - Math.Pow(eccentricity, 2)) - distance) / (eccentricity * distance));
            }
            else
            {
                angle = Math.Acos((semiMajorAxis * (1 - Math.Pow(eccentricity, 2)) - distance) / (eccentricity * distance));
            }

            float perpendicularVelocity = Convert.ToSingle(AngularMomentum / (Mass * distance));
            //float totalVelocitySquared = Convert.ToSingle((2 * TotalEnergy * distance + 2 * G * parentBody[2]) / (Mass * distance));
            float totalVelocitySquared = Convert.ToSingle(Math.Sqrt(G * parentBody[2] / distance));

            relativeVelocity = new Vector2(perpendicularVelocity, Convert.ToSingle(Math.Sqrt(Math.Pow(perpendicularVelocity, 2) - totalVelocitySquared)));

            angularVelocity = Math.Acos(AngularMomentum / Math.Sqrt(2 * TotalEnergy * Mass * Math.Pow(distance, 2) + 2 * G * parentBody[2] * Math.Pow(Mass, 2) * distance));

            radialPosition = new Vector2(Convert.ToSingle(distance), Convert.ToSingle(angle));

            radialPositionOut = radialPosition;

            convertToCartesian(parentBody, angle);
        }

        private void convertToCartesian(double[] parentBody, double angle)
        {
            if (angle < Math.PI / 2)
            {
                position.X = Convert.ToSingle(parentBody[0] - Math.Cos(radialPosition.Y) * radialPosition.X);
                position.Y = Convert.ToSingle(parentBody[1] - Math.Sin(radialPosition.Y) * radialPosition.X);
            }
            else if (angle < Math.PI)
            {
                position.X = Convert.ToSingle(parentBody[0] + Math.Cos(Math.PI / 2 - radialPosition.Y) * radialPosition.X);
                position.Y = Convert.ToSingle(parentBody[1] - Math.Sin(Math.PI / 2 - radialPosition.Y) * radialPosition.X);
            }
            else if (angle < Math.PI * 3 / 2)
            {
                position.X = Convert.ToSingle(parentBody[0] + Math.Cos(radialPosition.Y - Math.PI) * radialPosition.X);
                position.Y = Convert.ToSingle(parentBody[1] + Math.Sin(radialPosition.Y - Math.PI) * radialPosition.X);
            }
            else
            {
                position.X = Convert.ToSingle(parentBody[0] - Math.Cos(radialPosition.Y - Math.PI * 3 / 2) * radialPosition.X);
                position.Y = Convert.ToSingle(parentBody[1] + Math.Sin(radialPosition.Y - Math.PI * 3 / 2) * radialPosition.X);
            }
        }

        public Vector2 getPosition()
        {
            return position;
        }

        public Vector2 getVelocity()
        {
            return velocity;
        }

        public double getMass()
        {
            return Mass;
        }

        public double getTotalEnergy()
        {
            return TotalEnergy;
        }
    }
}
