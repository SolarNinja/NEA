using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Windows;
using System.Reflection.Metadata.Ecma335;
using System.DirectoryServices;
using System.Linq;
using System.Windows.Forms;
using System.Xml.Serialization;

namespace Orbital_Motion
{
    internal class Body
    {
        private int NullValue = Globals.NullValue;
        private Vector3 NullVector3 = Globals.NullVector3;

        private Vector3 relativePosition;
        private Vector3 relativeVelocity;
        private Vector3 relativePolarPosition;

        // Constants
        double Mass;
        double TotalEnergy, AngularMomentum;
        Vector3 SpecificAngularMomentum;
        static double G = 6.67 * Math.Pow(10, -11);

        // Psuedo-constants (constant unless a change in orbit occurs)
        //      - Longitude of Ascending Node and Periapsis has a reference direction of the relative positive x direction.
        private double eccentricity, semiMajorAxis, semiMinorAxis, periapsis, apoapsis, orbitalPeriod, timeOffset, longitudeOfAscendingNode, argumentOfPeriapsis, longitudeOfPeriapsis, inclination;
        private Dictionary<double, double[]> predeterminedPoints;


        /*   Properties of an orbiting body:
         *   Time: At what point in time is the body.
         *i  [Time offset: When the body was created; this needs to be subtracted from inputted time.
         *   Distance: How far from the orbiting body at a given time.
         * i [Start distance: How far it started from the parent body.
         *i       - Position can be turned into distance, angle from x and angle from y. 
         *   Angle: Angle from the true anomaly (given as relatively the positive x direction).
         * i [Angle offset: Given that the velocity does not start parallel to the tangent of the circle of its range and that if it does it is high enough, the body starts an an angle offset from the true anomaly.
         *        - Uncertain this needs to be explicit, as it may be implied by start speed and distance.
         *   Speed: How fast it is traveling at a given time.
         *   [Start speed: How fast it began.
         * i Speed's angle to the parallel: What direction the speed is pointing in relation to the tangent of the circle of its range.
         *i       - This can be found from the velocity vector, instead of using speed and angle.
         * i [Orbit angle offset: Angle the line connecting the apsides makes from the positive x axis, specifically the displacement vector from apoapsis to periapsis 
         * i [Longitude of ascending/decending node: The angle from a reference direction (presumably the positive x for the relative plane of the parent body) to the ascending/decending node. Decide on which to pick later.
         * i [Inclination: The angle the orbit makes from the horizontal plane.
         */

        // Unfortunately the one thing that cannot be used to create an ellipse is likely most useful; the orbital period. It's only equation is so complicated I do not what to attempt reversing it.
        // Currently without inclination from +x and xy plane.
        public Body(double[] parentBody, Vector3 startingRelativePosition, Vector3 startingRelativeVelocity,
                    double startingRelativeDistance, double startingRelativeSpeed, double longitudeOfAscendingNode,
                    double argumentOfPeriapsis, double longitudeOfPeriapsis, double inclination, double timeOffset,
                    double eccentricity, double semiMajorAxis, double semiMinorAxis, double periapsis, double apoapsis,
                    double Mass, double TotalEnergy, double startingKineticEnergy, double startingGravitationalPotentialEnergy,
                    double AngularMomentum, Vector3 SpecificAngularMomentum)
        {
            // This likely has redundant cases, so it may be worth double checking this later.
            // Also add checks later to make sure that all inputted data is consistent.
            void HandleInputs()
            {
                // A user may input various variables but together they do not uniquely define an orbit. This carries all possible variables that could be assigned too which, if inputted, would uniquely identify the orbit.
                string missingVariables = "";

                // This is the starting distance between the orbiting body and its parent. It starts null and is assigned to later.
                double startingRelativeAbsPosition = NullValue;

                // If either the starting position vector or the distance itself is given, the distance is assigned too accordingly.
                if (startingRelativePosition != NullVector3)
                {
                    startingRelativeAbsPosition = startingRelativePosition.Length();
                }
                else if (startingRelativeDistance != NullValue)
                {
                    startingRelativeAbsPosition = startingRelativeDistance;
                }

                // If the velocity exists, so does the speed. If not, it will be calculated later.
                double startingRelativeAbsVelocity = NullValue;
                if (startingRelativeVelocity != NullVector3)
                {
                    startingRelativeAbsVelocity = startingRelativeVelocity.Length();
                }
                else if (startingRelativeSpeed != NullValue)
                {
                    startingRelativeAbsVelocity = startingRelativeSpeed;
                }

                // Fix nulls later
                // Add error messages later and remember to clear errors between sections.
                // For things like TotalEnergy and AngularMomentum -1 is a possible value to equal and hence it may require a seperal NullValue.
                // Angular momentum doesn't include the right components of velocity, which should be fixed.

                // Orbital eccentricity is a unitless (dimensionless) parameter which gives a measure of how far an orbit is from a perfect circle.
                if (eccentricity != NullValue)
                {
                    this.eccentricity = eccentricity;
                }
                else
                {
                    // I ought to rewrite this horrendous if statement, it pains me. Not only that, it then repeats code.
                    if ((TotalEnergy != NullValue || ((startingGravitationalPotentialEnergy != NullValue || startingRelativeAbsPosition != NullValue) && (startingKineticEnergy != NullValue || startingRelativeAbsVelocity != NullValue))) && (AngularMomentum != NullValue || (startingRelativeAbsPosition != NullValue && startingRelativeAbsVelocity != NullValue)))
                    {
                        TotalEnergyAndAngularMomentum();
                        this.eccentricity = Math.Sqrt(1 + 2 * this.TotalEnergy * Math.Pow(this.AngularMomentum, 2) / (Math.Pow(Mass, 3) * Math.Pow(G * parentBody[2], 2))); ;
                    }
                    else if (semiMajorAxis != NullValue && semiMinorAxis != NullValue)
                    {
                        this.eccentricity = Math.Sqrt(1 - Math.Pow(semiMinorAxis / semiMajorAxis, 2));
                    }
                    else if (semiMajorAxis != NullValue)
                    {
                        if (periapsis != NullValue)
                        {
                            this.eccentricity = Math.Sqrt(1 - periapsis / semiMajorAxis);
                        }
                        else if (apoapsis != NullValue)
                        {
                            this.eccentricity = Math.Sqrt(apoapsis / semiMajorAxis - 1);
                        }
                    }
                    else
                    {
                        missingVariables += "";
                    }
                }

                // This is the distance from the centre of the orbit and the apoapsis or periapsis (the longest diameter).
                if (semiMajorAxis != NullValue)
                {
                    this.semiMajorAxis = semiMajorAxis;
                }
                else
                {
                    if (startingRelativeAbsPosition != NullValue && startingRelativeAbsVelocity != NullValue)
                    {
                        this.semiMajorAxis = 1 / (2 / startingRelativeAbsPosition - Math.Pow(startingRelativeAbsVelocity, 2) / (G * parentBody[2]));
                    }
                    else if (this.eccentricity != NullValue)
                    {
                        if (semiMinorAxis != NullValue)
                        {
                            this.semiMajorAxis = semiMinorAxis / Math.Sqrt(1 - Math.Pow(this.eccentricity, 2));
                        }
                        else if (periapsis != NullValue)
                        {
                            this.semiMajorAxis = periapsis / (1 - Math.Pow(this.eccentricity, 2));
                        }
                        else if (apoapsis != NullValue)
                        {
                            this.semiMajorAxis = apoapsis / (1 + Math.Pow(this.eccentricity, 2));
                        }
                        else
                        {
                            /*
                            missingVariables += "Semi-Minor Axis, Periapsis, Apoapsis, ";

                            if (startingRelativeDistance == NullValue && startingRelativeSpeed != NullValue)
                            {
                                missingVariables += "or Starting Velocity.";
                            }
                            else if (startingRelativeDistance != NullValue && startingRelativeSpeed == NullValue)
                            {
                                missingVariables += "or Starting Distance.";
                            }
                            else
                            {
                                missingVariables += "Starting Velocity, or Starting Distance.";
                            }

                            throw new Exception("You are missing either " + missingVariables);*/
                        }
                    }
                }

                // Given a semi-major axis, it is possible to define the second starting vector's magnitude.
                if (startingRelativePosition != NullVector3 && startingRelativeAbsVelocity == NullValue)
                {
                    startingRelativeAbsVelocity = Math.Sqrt(G * parentBody[2] * (2 / startingRelativeAbsPosition - 1 / this.semiMajorAxis));
                }
                else if (startingRelativeAbsPosition == NullValue && startingRelativeVelocity != NullVector3)
                {
                    startingRelativeAbsPosition = 2 * (1 / this.semiMajorAxis + Vector3.Dot(startingRelativeVelocity, startingRelativeVelocity));
                }

                // The semi-minor axis is the shortest diameter of an ellipse. It is perpendicular to the semi-major axis through the centre point if drawn on a graph.
                if (semiMinorAxis != NullValue)
                {
                    this.semiMinorAxis = semiMinorAxis;
                }
                else
                {
                    if (semiMajorAxis != NullValue && eccentricity != NullValue)
                    {
                        this.semiMinorAxis = this.semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity, 2));
                    }
                    else
                    {

                    }
                }

                // This is repeated (very unfortunately jankily) so I am putting it in a seperate nested function. I really should rewrite this.
                void TotalEnergyAndAngularMomentum()
                {
                    // Total energy is what it says on the tin. It will be conserved, and hence is vital to calulate position at any point.
                    if (TotalEnergy != NullValue)
                    {
                        this.TotalEnergy = TotalEnergy;
                    }
                    else
                    {
                        missingVariables += "Total Energy, ";

                        if (startingGravitationalPotentialEnergy != NullValue)
                        {
                            // Do nothing, it's fine.
                        }
                        else
                        {
                            missingVariables += "GPE, ";

                            if (startingRelativeAbsPosition != NullValue)
                            {
                                startingGravitationalPotentialEnergy = -1 * G * parentBody[2] * Mass / startingRelativeAbsPosition;
                            }
                            else
                            {
                                missingVariables += "or Starting Distance.";

                                throw new Exception("You are missing either " + missingVariables);
                            }
                        }

                        if (startingKineticEnergy != NullValue)
                        {
                            // Do nothing, it's fine.
                        }
                        else
                        {
                            missingVariables += "KE, ";

                            if (startingRelativeAbsVelocity != NullValue)
                            {
                                startingKineticEnergy = 0.5 * Mass * Vector3.Dot(startingRelativeVelocity, startingRelativeVelocity);
                            }
                            else
                            {
                                missingVariables += "or Starting Velocity.";

                                throw new Exception("You are missing either " + missingVariables);
                            }
                        }

                        this.TotalEnergy = startingKineticEnergy + startingGravitationalPotentialEnergy;
                        missingVariables = "";
                    }

                    // Angular momentum is the momentum in a given rotation. It is concerved, hence useful for calculations.
                    if (AngularMomentum != NullValue)
                    {
                        this.AngularMomentum = AngularMomentum;
                    }
                    else
                    {
                        missingVariables += "Angular Momentum, ";

                        if (startingRelativeAbsPosition != NullValue && startingRelativeAbsVelocity != NullValue)
                        {
                            this.AngularMomentum = Mass * startingRelativeAbsVelocity * startingRelativeAbsPosition;
                            missingVariables = "";
                        }
                        else
                        {
                            if (startingRelativeAbsPosition == NullValue && startingRelativeAbsVelocity != NullValue)
                            {
                                missingVariables += "or Starting Velocity.";
                            }
                            else if (startingRelativeAbsPosition != NullValue && startingRelativeAbsVelocity == NullValue)
                            {
                                missingVariables += "or Starting Distance.";
                            }
                            else
                            {
                                missingVariables += "Starting Velocity, or Starting Distance.";
                            }

                            throw new Exception("You are missing either " + missingVariables);
                        }
                    }
                }

                TotalEnergyAndAngularMomentum();

                // This is the shortest distance from the parent body.
                if (periapsis != NullValue)
                {
                    this.periapsis = periapsis;
                }
                else
                {
                    this.periapsis = this.semiMajorAxis * (1 - this.eccentricity) + 0.01;
                }

                // This is the largest distance from the parent body.
                if (apoapsis != NullValue)
                {
                    this.apoapsis = apoapsis;
                }
                else
                {
                    this.apoapsis = this.semiMajorAxis * (1 + this.eccentricity) - 0.01;
                }

                // This is the vector perpendicular from the other 2 vectors which is a measure of angular momentum / mass.
                if (SpecificAngularMomentum != NullVector3)
                {
                    this.SpecificAngularMomentum = SpecificAngularMomentum;
                }
                else if (startingRelativePosition != NullVector3 && startingRelativeVelocity != NullVector3)
                {
                    this.SpecificAngularMomentum = Vector3.Cross(startingRelativePosition, startingRelativeVelocity);
                }
                else if (this.AngularMomentum != NullValue)
                {
                    this.SpecificAngularMomentum = new Vector3(0, 0, Convert.ToSingle(this.AngularMomentum / Mass));
;               }

                // The angle between the equatorial plane (in this case xy) and the plane of the orbit.
                if (inclination != NullValue)
                {
                    this.inclination = inclination;
                }
                else
                {
                    this.inclination = Math.Acos(this.SpecificAngularMomentum.Z / this.SpecificAngularMomentum.Length());
                }

                // n is the vector pointing towards the ascending node. If the ascending node doesn't exist (inclination = 0 or pi) it is assigned as the nullVector3. It is used for the following equations.
                Vector3 n = NullVector3;
                if (inclination != 0 && inclination != Math.PI)
                {
                    n = Vector3.Cross(new Vector3(0, 0, 1), this.SpecificAngularMomentum);
                }

                // e is the eccentricity vector, a vector with no unit that points from apoapsis to periapsis with magnitude equal to scalar eccentricity.
                // It is important to note this can only be calculated knowing the starting vectors, but this can still be bipassed by knowing the LoP/ 
                // This is also used for the following equations.
                Vector3 e = NullVector3;
                if (startingRelativePosition != NullVector3 && startingRelativeVelocity != NullVector3)
                {
                    e = Vector3.Divide(Vector3.Cross(startingRelativeVelocity, this.SpecificAngularMomentum), Convert.ToSingle(G * parentBody[2])) - Vector3.Divide(startingRelativePosition, Convert.ToSingle(startingRelativeDistance));
                }

                // The angle between the +x direction and the ascending node.
                if (longitudeOfAscendingNode != NullValue)
                {
                    this.longitudeOfAscendingNode = longitudeOfAscendingNode;
                }
                else
                {
                    if (inclination == 0)
                    {
                        this.longitudeOfAscendingNode = 0;
                    }
                    else
                    {
                        if (n != nullVector3)
                        if (n.Y >= 0)
                        {
                            this.longitudeOfAscendingNode = Math.Acos(n.X / n.Length());
                        }
                        else
                        {
                            this.longitudeOfAscendingNode = 2 * Math.PI - Math.Acos(n.X / n.Length());
                        }
                    }
                }

                // The angle between the +x direction and the periapsis. This is being used as if inclination == 0 the argument of periapsis defaults to 0 w
                if (argumentOfPeriapsis != NullValue)
                {
                    this.argumentOfPeriapsis = argumentOfPeriapsis;
                }
                else if (startingRelativePosition != NullVector3 && startingRelativeVelocity != NullVector3)
                {
                    if (inclination == 0)
                    {
                        this.argumentOfPeriapsis = 0;
                    }
                    else if (longitudeOfPeriapsis != NullValue && this.longitudeOfAscendingNode != NullValue)
                    {
                        this.argumentOfPeriapsis = longitudeOfPeriapsis - this.longitudeOfAscendingNode;
                    }
                    else
                    {
                        this.argumentOfPeriapsis = Math.Acos(Vector3.Dot(n, e) / (n.Length() * e.Length()));

                        if (e.Z < 0)
                        {
                            this.argumentOfPeriapsis = 2 * Math.PI - this.argumentOfPeriapsis;
                        }
                    }
                }

                if (longitudeOfPeriapsis != NullValue)
                {
                    this.longitudeOfPeriapsis = longitudeOfPeriapsis;
                }
                else 
                {
                    this.longitudeOfPeriapsis = this.longitudeOfAscendingNode + this.argumentOfPeriapsis; 
                }

                relativePosition = startingRelativePosition;
                relativeVelocity = startingRelativeVelocity;
                this.timeOffset = timeOffset;
                this.Mass = Mass;
            }

            HandleInputs();

            //predeterminedPoints = createPoints(10000, parentBody, periapsis, apoapsis, 0, out orbitalPeriod);
        }

        /*public Body(Vector2 startingRelativePosition, Vector2 startingVelocity, double Mass, double[] parentBody)
        {
            // Are only used here so are only defined locally.
            double startingKineticEnergy, startingGravitationalPotentialEnergy;

            //Assigning starting values from the function inputs.
            relativePosition = startingRelativePosition;
            relativeVelocity = startingVelocity;
            this.Mass = Mass;

            // A value of r to use to gather the psuedo-constant's values.
            startingDistance = Math.Sqrt(Math.Pow(relativePosition.X, 2) + Math.Pow(relativePosition.Y, 2));

            // Energies at start and the total; please note the use of TE = 1/2 * GPE only works for a perfect circular orbit, hence cannot be used.
            startingGravitationalPotentialEnergy = -1 * G * parentBody[2] * Mass / startingDistance;
            startingKineticEnergy = 0.5 * Mass * Vector2.Dot(relativeVelocity, relativeVelocity);
            TotalEnergy = startingGravitationalPotentialEnergy + startingKineticEnergy;
            
            // Defining another psuedo-constant.
            AngularMomentum = Mass * relativeVelocity.Length() * startingDistance;

            eccentricity = Math.Sqrt(1 + 2 * TotalEnergy * Math.Pow(AngularMomentum, 2) / (Math.Pow(Mass, 3) * Math.Pow(G * parentBody[2], 2)));
            semiMajorAxis = 1 / (2 / startingDistance - Vector2.Dot(relativeVelocity,relativeVelocity) / (G * parentBody[2]));
            semiMinorAxis = semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity, 2));
            // These are minorly adjusted to make sure no errors occur when calculating the extreme cases.
            periapsis = semiMajorAxis * (1 - eccentricity) + 0.01;
            apoapsis = semiMajorAxis * (1 + eccentricity) - 0.01;

            // Orbital period is defined here due to being more convienient in timing with the maths.
            predeterminedPoints = createPoints(10000, parentBody, periapsis, apoapsis, 0, out orbitalPeriod);
        }*/

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

        /*public void returnValues(double time, double[] parentBody, out Vector2 radialPositionOut, out Vector2 relativeVelocity, out double angularVelocity)
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

            //relativePolarPosition = new Vector2(Convert.ToSingle(distance), Convert.ToSingle(angle));

            //radialPositionOut = relativePolarPosition;

            convertToCartesian(parentBody, angle);
        }*/

        private void convertToCartesian(double[] parentBody, double angle)
        {
            relativePosition.X = Convert.ToSingle(parentBody[0] - Math.Cos(relativePolarPosition.Y) * relativePolarPosition.X);
            relativePosition.Y = Convert.ToSingle(parentBody[1] - Math.Sin(relativePolarPosition.Y) * relativePolarPosition.X);
        }
    }
}


/*         
            this.relativePosition = relativePosition;
            this.velocity = velocity;
            this.radialPosition = radialPosition;
            this.startingDistance = startingDistance;
            Mass = mass;
            this.eccentricity = eccentricity;
            this.semiMajorAxis = semiMajorAxis;
            this.semiMinorAxis = semiMinorAxis;
            this.periapsis = periapsis;
            this.apoapsis = apoapsis;
            this.orbitalPeriod = orbitalPeriod;
            this.predeterminedPoints = predeterminedPoints;
*/