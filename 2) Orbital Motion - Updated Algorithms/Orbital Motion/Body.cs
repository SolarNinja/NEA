using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Windows;
using System.Reflection.Metadata.Ecma335;
using System.DirectoryServices;
using System.Linq;

namespace Orbital_Motion
{
    internal class Body
    {
        private Vector2 position;
        private Vector2 velocity;
        private Vector2 radialPosition;

        private double kineticEnergy, gravitationalPotentialEnergy, startingDistance;

        // Constants
        double Mass;
        static double TotalEnergy, AngularMomentum;
        //static double G = 6.6743 * Math.Pow(10, -11);
        static double G = 6.67 * Math.Pow(10, -11);

        // Psuedo constants (constant unless a change in orbit occurs)
        private double eccentricity, semiMajorAxis, semiMinorAxis, periapsis, apoapsis;
        private Dictionary<double, double[]> predeterminedPoints;

        public Body(Vector2 position, Vector2 velocity, double Mass, double[] parentBody)
        {
            this.position = position;
            this.velocity = velocity;
            this.Mass = Mass;

            startingDistance = Math.Sqrt(Math.Pow(position.X - parentBody[0], 2) + Math.Pow(position.Y - parentBody[1], 2));

            gravitationalPotentialEnergy = -1 * G * parentBody[2] * Mass / startingDistance;
            kineticEnergy = 0.5 * Mass * Vector2.Dot(velocity, velocity);

            TotalEnergy = gravitationalPotentialEnergy + kineticEnergy;
            AngularMomentum = Mass * velocity.Length() * startingDistance;

            eccentricity = Math.Sqrt(1 + 2 * TotalEnergy * Math.Pow(AngularMomentum, 2) / (Math.Pow(Mass, 3) * Math.Pow(G * parentBody[2], 2)));
            semiMajorAxis = 1 / (2 / startingDistance - Vector2.Dot(velocity,velocity) / (G * parentBody[2]));
            semiMinorAxis = semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity, 2));
            // These are minorly adjusted to make sure no errors occur when calculating the extreme cases.
            periapsis = semiMajorAxis * (1 - eccentricity) + 0.01;
            apoapsis = semiMajorAxis * (1 + eccentricity) - 0.01;

            predeterminedPoints = createPoints(10000, parentBody, periapsis, apoapsis);
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
        private Dictionary<double, double[]> createPoints(int numberOfIntervals, double[] parentBody, double periapsis, double apoapsis)
        {
            Dictionary<double, double[]> predeterminedPoints = new Dictionary<double, double[]>();

            double range = (apoapsis - periapsis) / 2;

            double V = TotalEnergy / Mass;
            double A = AngularMomentum / Mass;
            double K = Math.Pow(G * parentBody[2], 2) + 2 * V * Math.Pow(A, 2);

            for (int i = 0; i < numberOfIntervals; i++)
            {
                double R = periapsis + i * range / numberOfIntervals;
                double U = 2 * V * Math.Pow(R, 2) + 2 * G * parentBody[2] * R - Math.Pow(A, 2);
                double time = -1 * G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Sqrt(-2 * V * U / K))) + 1 / (2 * V) * Math.Sqrt(U);
                predeterminedPoints.Add(i, new double[] {R, time});
            }

            return predeterminedPoints;
        }

        public void returnValues(double time, double[] parentBody, out Vector2 radialPositionOut, out Vector2 relativeVelocity, out double angularVelocity)
        {
            /*int count = 0;
            Complex func = 0;
            void newtonRaphsonDistance(double time, Complex startingR, out double distance)
            {
                Complex R = startingR;

                Complex U = 2 * V * Complex.Pow(R, 2) + 2 * G * parentBody[2] * R - Math.Pow(A, 2);
                // double i = Math.Sqrt(-2 * V * U / K);
                func = -1 * G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Complex.Asin(Complex.Sqrt(-2 * V * U / K))) + 1 / (2 * V) * Complex.Sqrt(U) - time;
                // G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Sqrt(-2 * V * U / K))) + 1 / (2 * V) * Math.Sqrt(U) - time;
                Complex funcDerivative = R / Complex.Sqrt(U);

                Complex newtonRaphsonAnswer = R - func / funcDerivative;

                if (count < 100)
                {
                    count++;
                    newtonRaphsonDistance(time, newtonRaphsonAnswer, out distance);
                }
                distance = newtonRaphsonAnswer.Magnitude;
            }

            void fixedPointIteration(double time, double startingR, out double distance)
            {
                double R = startingR;

                double U = 2 * V * Math.Pow(R, 2) + 2 * G * parentBody[2] * R - Math.Pow(A, 2);
                //func = -1 * G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Sqrt(-2 * V * U / K))) + 1 / (2 * V) * Math.Sqrt(U) - time;
                double iteration = (Math.Pow(2 * V * time - G * parentBody[2] / 2 * (Math.Sqrt(-2 / V) * Math.Asin(Math.Sqrt(-2 * V * U / K))), 2) - 2 * V * Math.Pow(R, 2) + Math.Pow(A, 2)) / (2 * G * parentBody[2]);
                
                if (count < 100)
                {
                    count++;
                    fixedPointIteration(time, iteration, out distance);
                }
                distance = iteration;
            }*/

            double predeterminedMethod()
            {
                double[] predeterminedAdjustedTimes = new double[predeterminedPoints.Count];

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

            double angle = Math.Acos((semiMajorAxis * (1 - Math.Pow(eccentricity, 2)) - distance) / (eccentricity * distance));

            float perpendicularVelocity = Convert.ToSingle(AngularMomentum / (Mass * distance));
            float totalVelocitySquared = Convert.ToSingle((2 * TotalEnergy * distance + 2 * G * parentBody[2]) / (Mass * distance));

            relativeVelocity = new Vector2(perpendicularVelocity, Convert.ToSingle(Math.Sqrt(Math.Pow(perpendicularVelocity, 2) - totalVelocitySquared)));

            angularVelocity = Math.Acos(AngularMomentum / Math.Sqrt(2 * TotalEnergy * Mass * Math.Pow(distance, 2) + 2 * G * parentBody[2] * Math.Pow(Mass, 2) * distance));

            radialPosition = new Vector2(Convert.ToSingle(distance), Convert.ToSingle(angle));

            radialPositionOut = radialPosition;

            convertToCartesian(parentBody);
        }

        private void convertToCartesian(double[] parentBody)
        {
            position.X = Convert.ToSingle(parentBody[0] - Math.Cos(radialPosition.Y) * radialPosition.X);
            position.Y = Convert.ToSingle(parentBody[1] - Math.Sin(radialPosition.Y) * radialPosition.X);
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
