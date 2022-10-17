using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Windows;
using System.Reflection.Metadata.Ecma335;

namespace Orbital_Motion
{
    internal class Body
    {
        private Vector2 position;
        private Vector2 velocity;
        private double mass, kineticEnergy, gravitationalPotentialEnergy, startingDistance;

        // Constants
        static double TotalEnergy, AngularMomentum;
        static double G = 6.6743 * Math.Pow(10, -11);

        public Body(Vector2 position, Vector2 velocity, double mass, double[] parentBody)
        {
            this.position = position;
            this.velocity = velocity;
            this.mass = mass;

            startingDistance = Math.Sqrt(Math.Pow(position.X - parentBody[0], 2) + Math.Pow(position.Y - parentBody[1], 2));

            gravitationalPotentialEnergy = -1 * G * parentBody[2] * mass / startingDistance;
            kineticEnergy = 0.5 * Vector2.Dot(velocity, velocity);

            TotalEnergy = gravitationalPotentialEnergy + kineticEnergy;
            AngularMomentum = mass * velocity.Length() * startingDistance;
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
            double l = mass * (velocity.X * position.Y - velocity.Y * position.X);

            // The eccentricity.
            double e = Math.Sqrt(1 + (2 * TotalEnergy * l * l) / (mass * mass * mass * mew * mew));

            // Interestingly, the eccentricity vector of a craft of position p and velocity v is given by a different equation.
            Vector2 eVector = Vector2.Divide(velocity * velocity * position, Convert.ToSingle(mew)) - Vector2.Divide(position, Convert.ToSingle(position.Length()));

            // The true anomaly, as calculated using the eccentricity vector.
            double theta = Math.Acos(Vector2.Dot(eVector, position) / (eVector.Length() * position.Length()));

            // The range between the body and the parent.
            double r = (l * l) / (mass * mass * mew) * 1 / (1 + e * Math.Cos(theta));

            return r;
        }

        /*public void returnValues(int time, out int angle, out int distance, out Vector2 relativeVelocity, out int angularVelocity)
        {

        }*/

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
            return mass;
        }

        public double getTotalEnergy()
        {
            return TotalEnergy;
        }
    }
}
