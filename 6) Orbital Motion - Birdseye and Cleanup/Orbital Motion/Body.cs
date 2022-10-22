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
using System.Drawing.Drawing2D;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TaskbarClock;
using System.Drawing;

namespace Orbital_Motion
{
    internal class Body
    {
        // This code is a right mess.

        record pointInformation(double time, double angularVelocity, double trueAnomaly, Vector3 relativePosition, Vector3 relativeVelocity);

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
        private pointInformation[] predeterminedPoints;

        // Unfortunately the one thing that cannot be used to create an ellipse is likely most useful; the orbital period. It's only equation is so complicated I do not what to attempt reversing it.
        // Currently without inclination from +x and xy plane.
        public Body(double[] parentBody, double Mass, double timeOffset, Vector3 startingRelativePosition, 
                    Vector3 startingRelativeVelocity, double startingRelativeDistance, double startingRelativeSpeed,
                    double startingAngleOfSpeed, double longitudeOfAscendingNode, double argumentOfPeriapsis, double longitudeOfPeriapsis, 
                    double inclination, double eccentricity, double semiMajorAxis, double semiMinorAxis, 
                    double periapsis, double apoapsis, double TotalEnergy, double startingKineticEnergy,
                    double startingGravitationalPotentialEnergy, double AngularMomentum, Vector3 SpecificAngularMomentum)
        {
            // This likely has redundant cases, so it may be worth double checking this later.
            // Also add checks later to make sure that all inputted data is consistent.
            double startingPerpendicularAbsVelocity = NullValue;
            double startingRelativeAbsPosition = NullValue;
            void HandleInputs()
            {
                // A user may input various variables but together they do not uniquely define an orbit. This carries all possible variables that could be assigned too which, if inputted, would uniquely identify the orbit.
                string missingVariables = "";

                // This is the starting distance between the orbiting body and its parent. It starts null and is assigned to later.
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

                // This is the vector perpendicular from the other 2 vectors which is a measure of angular momentum / mass.
                if (SpecificAngularMomentum != NullVector3)
                {
                    this.SpecificAngularMomentum = SpecificAngularMomentum;
                }
                else if (startingRelativePosition != NullVector3 && startingRelativeVelocity != NullVector3)
                {
                    this.SpecificAngularMomentum = Vector3.Cross(startingRelativePosition, startingRelativeVelocity);
                }

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
                if (this.inclination != 0 && this.inclination != Math.PI)
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
                    if (this.inclination == 0 || this.inclination == Math.PI)
                    {
                        this.longitudeOfAscendingNode = 0;
                    }
                    else
                    {
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
                else 
                {
                    if (this.inclination == 0 || this.inclination == Math.PI)
                    {
                        // Given an inclination of 0 this will be reassigned to later, given its use of eccentricity.
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

                // Maybe make angular momentum reversable as well?
                if (startingAngleOfSpeed != NullValue)
                {
                    startingPerpendicularAbsVelocity = startingRelativeAbsVelocity * Math.Cos(startingAngleOfSpeed);
                }
                else if (startingRelativeVelocity != NullVector3 && startingRelativePosition != NullVector3)
                {
                    if (this.inclination != 0)
                    {
                        // This method quite possibly is overcomplicating the maths. I'll need to revist this later.
                        Vector3 flattenedRelativeStartPosition = startingRelativePosition;
                        Vector3 flattenedRelativeStartVelocity = startingRelativeVelocity;

                        // This rotates the line between the ascending and decending node to be along the positive x axis.
                        double angleFromAscendingNode = 2 * Math.PI - this.longitudeOfAscendingNode;
                        Matrix3x3 rotationMatrixForAscendingNode = new Matrix3x3(new double[,] { { Math.Cos(angleFromAscendingNode), Math.Sin(angleFromAscendingNode), 0 }, { -Math.Sin(angleFromAscendingNode), Math.Cos(angleFromAscendingNode), 0 }, { 0, 0, 1 } });

                        flattenedRelativeStartPosition = rotationMatrixForAscendingNode.Transform(flattenedRelativeStartPosition);
                        flattenedRelativeStartVelocity = rotationMatrixForAscendingNode.Transform(flattenedRelativeStartVelocity);

                        // This rotates the orbit along it's node line (the x axis), hence giving a flat orbit.
                        double angleFromInclination = 2 * Math.PI - this.inclination;
                        Matrix3x3 rotationMatrixForInclination = new Matrix3x3(new double[,] { { 1, 0, 0 }, { 0, Math.Cos(angleFromInclination), Math.Sin(angleFromInclination) }, { 0, -Math.Sin(angleFromInclination), Math.Cos(angleFromInclination) } });
                        
                        flattenedRelativeStartPosition = rotationMatrixForInclination.Transform(flattenedRelativeStartPosition);
                        flattenedRelativeStartVelocity = rotationMatrixForInclination.Transform(flattenedRelativeStartVelocity);

                        // This rotates the orbit such that its periapsis and apoapsis lie on the x axis, specifically such that the periapsis is furthest along the positive x direction.
                        double angleFromPeriapsis = 2 * Math.PI - this.argumentOfPeriapsis;
                        Matrix3x3 rotationMatrixForPeriapsis = new Matrix3x3(new double[,] { { Math.Cos(angleFromPeriapsis), Math.Sin(angleFromPeriapsis), 0 }, { -Math.Sin(angleFromPeriapsis), Math.Cos(angleFromPeriapsis), 0 }, { 0, 0, 1 } });

                        flattenedRelativeStartPosition = rotationMatrixForPeriapsis.Transform(flattenedRelativeStartPosition);
                        flattenedRelativeStartVelocity = rotationMatrixForPeriapsis.Transform(flattenedRelativeStartVelocity);

                        startingAngleOfSpeed = Math.PI / 2 + Math.Atan(flattenedRelativeStartPosition.Y / flattenedRelativeStartPosition.X) - Math.Atan(flattenedRelativeStartVelocity.Y / flattenedRelativeStartVelocity.X);

                        startingPerpendicularAbsVelocity = startingRelativeAbsVelocity * Math.Cos(startingAngleOfSpeed);
                    }
                    else
                    {
                        startingAngleOfSpeed = Math.PI / 2 + Math.Atan(startingRelativePosition.Y / startingRelativePosition.X) - Math.Atan(startingRelativeVelocity.Y / startingRelativeVelocity.X);

                        startingPerpendicularAbsVelocity = startingRelativeAbsVelocity * Math.Cos(startingAngleOfSpeed);
                    }
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
                            if (startingRelativeAbsPosition == NullValue)
                            {
                                startingRelativeAbsPosition = -1 * G * parentBody[2] * Mass / startingGravitationalPotentialEnergy;
                            }
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
                            if (startingRelativeAbsVelocity == NullValue)
                            {
                                startingRelativeAbsVelocity = Math.Sqrt(2 * startingKineticEnergy / Mass);
                            }
                        }
                        else
                        {
                            missingVariables += "KE, ";

                            if (startingRelativeAbsVelocity != NullValue)
                            {
                                startingKineticEnergy = 0.5 * Mass * Math.Pow(startingRelativeAbsVelocity, 2);
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

                    // Angular momentum is the momentum in a given rotation. It is conserved, hence useful for calculations.
                    if (AngularMomentum != NullValue)
                    {
                        this.AngularMomentum = AngularMomentum;

                        if (startingRelativeAbsVelocity == NullValue && startingPerpendicularAbsVelocity == NullValue && startingAngleOfSpeed != NullValue)
                        {
                            startingPerpendicularAbsVelocity = AngularMomentum / (Mass * startingRelativeAbsPosition);
                            startingRelativeAbsVelocity = startingPerpendicularAbsVelocity / Math.Cos(startingAngleOfSpeed);
                        }
                    }
                    else
                    {
                        missingVariables += "Angular Momentum, ";

                        if (startingRelativeAbsPosition != NullValue && startingRelativeAbsVelocity != NullValue)
                        {
                            this.AngularMomentum = Mass * startingPerpendicularAbsVelocity * startingRelativeAbsPosition;
                            missingVariables = "";
                        }
                        else
                        {
                            if (startingRelativeAbsPosition == NullValue && startingPerpendicularAbsVelocity != NullValue)
                            {
                                missingVariables += "or Starting Velocity.";
                            }
                            else if (startingRelativeAbsPosition != NullValue && startingPerpendicularAbsVelocity == NullValue)
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

                relativePosition = startingRelativePosition;
                relativeVelocity = startingRelativeVelocity;
                this.timeOffset = timeOffset;
                this.Mass = Mass;
            }

            HandleInputs();

            predeterminedPoints = createPoints(10000, parentBody, startingRelativeAbsPosition, ref orbitalPeriod);
        }

        private pointInformation[] createPoints(int numberOfIntervals, double[] parentBody, double startDistance, ref double orbitalPeriod)
        {
            Dictionary<int, double[]> predeterminedDistances = new Dictionary<int, double[]>();

            // V is the specific energy.
            double V = TotalEnergy / Mass;

            // A is the absolute value of the specific angular momentum.
            double A = AngularMomentum / Mass;

            // K is a constant used in the following equations.
            double K = Math.Pow(G * parentBody[2], 2) + 2 * V * Math.Pow(A, 2);

            // midpointDistance is the distance halfway between apoapsis and periapsis.
            double midpointDistance = -G * parentBody[2] / (V * 2) + 0.01;

            // U is a variable determined by r; midpointAdjustedU is hence a variable determined by midpointDistance.
            double midpointAdjustedU = 2 * V * Math.Pow(midpointDistance, 2) + 2 * G * parentBody[2] * midpointDistance - Math.Pow(A, 2);

            // U is a variable determined by r; periapsisU is hence a variable determined by periapsis.
            double periapsisU = 2 * V * Math.Pow(periapsis, 2) + 2 * G * parentBody[2] * periapsis - Math.Pow(A, 2);

            // midpointTimeAdjustment is the time added to the second equation for time too adjust for the discontinuity in the graphs.
            double midpointTimeAdjustment = -2 * G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * midpointAdjustedU / K), 5)));

            // orbitalPeriod is the total time taken to complete one orbit.
            orbitalPeriod = 2 * (G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * periapsisU / K), 5))) + 1 / (2 * V) * Math.Sqrt(periapsisU) + midpointTimeAdjustment);

            double range = apoapsis - periapsis;

            // Given that the position the body starts at is not the apoapsis, there is an offset in the time from the natural time given from starting at apoapsis.
            double startPositionTimeOffset = 0;
            // When calculating the offset at apoapsis or periapsis please slightly align based on them to not cause errors? + for peri - for apo
            CalculatePoint(startDistance - 0.0001, startPositionTimeOffset, orbitalPeriod, out startPositionTimeOffset);

            // When at a non-inclined orbit, standard methods of calculating the argument of periapsis do not work due to their reliance on the ascending node's position.
            // A different method can be used; given a starting position you can calculate it's angle from the positive x direction as well as its true anomaly, hence finding the difference to calculate the overall offset.

            // Check whether this deals with +-
            double startingTrueAnomaly = Math.Acos((semiMajorAxis * (1 - Math.Pow(eccentricity, 2)) - startDistance) / (eccentricity * startDistance)); ;
            if (inclination == 0 || inclination == 2 * Math.PI)
            {
                double acuteAngle = Math.Atan(Math.Abs(relativePosition.Y / relativePosition.X));
                if (relativePosition.X >= 0 && relativePosition.Y >= 0)
                {
                    argumentOfPeriapsis = acuteAngle - startingTrueAnomaly;
                }
                else if (relativePosition.X < 0 && relativePosition.Y >= 0)
                {
                    argumentOfPeriapsis = Math.PI - acuteAngle - startingTrueAnomaly;
                }
                else if (relativePosition.X < 0 && relativePosition.Y < 0)
                {
                    argumentOfPeriapsis = Math.PI + acuteAngle - startingTrueAnomaly;
                }
                else if (relativePosition.X >= 0 && relativePosition.Y < 0)
                {
                    argumentOfPeriapsis = 2 * Math.PI - acuteAngle - startingTrueAnomaly;
                }
            }

            // Might work at half cant remember test later.
            for (int i = 0; i < 2 * numberOfIntervals; i += 2)
            {
                double R = periapsis + i * range / (2 * numberOfIntervals);
                double time;
                CalculatePoint(R, startPositionTimeOffset, orbitalPeriod, out time);

                predeterminedDistances.Add(i, new double[] { R, time });
                predeterminedDistances.Add(i + 1, new double[] { R, orbitalPeriod - time });
            }

            pointInformation[] predeterminedPoints = new pointInformation[predeterminedDistances.Count];
            foreach (var point in predeterminedDistances)
            {
                double time = point.Value[1];

                double trueAnomaly = Math.Acos((semiMajorAxis * (1 - Math.Pow(eccentricity, 2)) - point.Value[0]) / (eccentricity * point.Value[0]));

                double angularVelocity = Math.Acos(AngularMomentum / Math.Sqrt(2 * TotalEnergy * Mass * Math.Pow(point.Value[0], 2) + 2 * G * parentBody[2] * Math.Pow(Mass, 2) * point.Value[0]));

                Vector3 relativePosition = new Vector3(Convert.ToSingle(point.Value[0] * Math.Cos(trueAnomaly)), Convert.ToSingle(point.Value[0] * Math.Sin(trueAnomaly)), 0);

                float perpendicularVelocity = Convert.ToSingle(AngularMomentum / (Mass * point.Value[0]));
                float totalVelocitySquared = Convert.ToSingle((2 * TotalEnergy * point.Value[0] + 2 * G * parentBody[2]) / (Mass * point.Value[0]));
                //float totalVelocitySquared = Convert.ToSingle(Math.Sqrt(G * parentBody[2] / pointInformation[0]));
                Vector3 relativeVelocity = new Vector3(perpendicularVelocity, Convert.ToSingle(Math.Sqrt(Math.Pow(perpendicularVelocity, 2) - totalVelocitySquared)), 0);

                RotatePoint(ref relativePosition, ref relativeVelocity);

                predeterminedPoints[point.Key] = new pointInformation(time, angularVelocity, trueAnomaly, relativePosition, relativeVelocity);
            }

            return predeterminedPoints;

            void CalculatePoint(double R, double startPositionTimeOffset, double orbitalPeriod, out double time)
            {
                double U = 2 * V * Math.Pow(R, 2) + 2 * G * parentBody[2] * R - Math.Pow(A, 2);
                time = 0;

                // I actually believe this to be more efficient than defining repeating parts as it saves on read-writes.
                // Also note the use of the round; this is because an error would occur at halfway where, due to precision errors, the arcsin would recieve a value ever so slightly above 1. I might need to revisit this error to make it's handling more efficient.
                if (R < midpointDistance)
                {
                    time = Math.Abs(-1 * G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + timeOffset - startPositionTimeOffset) % orbitalPeriod;
                }
                else if (R >= midpointDistance)
                {
                    time = Math.Abs(G * parentBody[2] / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + midpointTimeAdjustment + timeOffset - startPositionTimeOffset) % orbitalPeriod;
                }
            }

            void RotatePoint(ref Vector3 relativePosition, ref Vector3 relativeVelocity)
            {
                // This rotates the orbit such that its periapsis is at the correct angle to the positive x direction.
                Matrix3x3 rotationMatrixForPeriapsis = new Matrix3x3(new double[,] { { Math.Cos(argumentOfPeriapsis), Math.Sin(argumentOfPeriapsis), 0 }, { -Math.Sin(argumentOfPeriapsis), Math.Cos(argumentOfPeriapsis), 0 }, { 0, 0, 1 } });

                relativePosition = rotationMatrixForPeriapsis.Transform(relativePosition);
                relativeVelocity = rotationMatrixForPeriapsis.Transform(relativeVelocity);

                // This rotates the orbit to be the correct inclination.
                Matrix3x3 rotationMatrixForInclination = new Matrix3x3(new double[,] { { 1, 0, 0 }, { 0, Math.Cos(inclination), Math.Sin(inclination) }, { 0, -Math.Sin(inclination), Math.Cos(inclination) } });

                relativePosition = rotationMatrixForInclination.Transform(relativePosition);
                relativeVelocity = rotationMatrixForInclination.Transform(relativeVelocity);

                // This rotates the orbit such that the ascending node is at the correct angle to the positive x direction.
                Matrix3x3 rotationMatrixForAscendingNode = new Matrix3x3(new double[,] { { Math.Cos(longitudeOfAscendingNode), Math.Sin(longitudeOfAscendingNode), 0 }, { -Math.Sin(longitudeOfAscendingNode), Math.Cos(longitudeOfAscendingNode), 0 }, { 0, 0, 1 } });

                relativePosition = rotationMatrixForAscendingNode.Transform(relativePosition);
                relativeVelocity = rotationMatrixForAscendingNode.Transform(relativeVelocity);
            }
        }

        // When calculating position, rotate based on argument of periapsis then inclination then longitude of ascending node
        // For an equatorial or retrograde orbit, calculate longitude of periapsis through a different method and rotate once based on that.
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