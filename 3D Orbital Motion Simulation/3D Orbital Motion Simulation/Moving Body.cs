using System.Drawing;
using System.Numerics;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.Xml;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TaskbarClock;

namespace _3D_Orbital_Motion_Simulation
{
    internal class MovingBody : Body
    {
        // **The nullability of these are likely unneccessary due to the nature of creating a body requiring automated processes and not nice code writing.
        public record PointInformation(double time = double.NaN, double angleOfVelocityToTangent = double.NaN, double trueAnomaly = double.NaN, Vector3 relativePosition = new Vector3(), Vector3 relativeVelocity = new Vector3());
        public record OrbitInformation(double longitudeOfAscendingNode, double argumentOfPeriapsis, double longitudeOfPeriapsis, double inclination, double eccentricity, double semiMajorAxis, double semiMinorAxis, double periapsis, double apoapsis, double AngularMomentum, Vector3 SpecificAngularMomentum, Vector3 eccentricityVector, double TotalEnergy, double hillSphereRadius, double orbitalPeriod);

        // The gravitational constant; used for maths. **May be moveable into the orbit creation function.
        static double G = 6.67 * Math.Pow(10, -11);

        // The body's name.
        public override string name { get; set; }

        // The body's mass.
        public override double Mass { get; set; }

        private int startingTimeFromEpoch { get; set; }

        // The information about a body's orbit.
        public OrbitInformation orbitInformation { get; private set; }

        // The body's current point information.
        public PointInformation currentPoint { get; private set; }

        // The full orbit's points precalculated on initialisation.
        private PointInformation[] predeterminedPoints;

        public MovingBody(ref bool failed, Tree tree, string name, double Mass, double ParentMass, double parentHillSphereRadius, int startingTimeFromEpoch, Vector3 startingRelativePosition, Vector3 startingRelativeVelocity, Vector3 startingAbsolutePosition)
        { 
            this.name = name;

            this.Mass = Mass;

            this.startingTimeFromEpoch = startingTimeFromEpoch;

            failed = InitialiseOrbitConstantsFromStateVectors(tree, ParentMass, parentHillSphereRadius, startingRelativePosition, startingRelativeVelocity, startingAbsolutePosition);
        }

        public MovingBody(ref bool failed, Tree tree, string name, double Mass, double ParentMass, double parentHillSphereRadius, Vector3 parentPosition, int startingTimeFromEpoch, double startingTrueAnomaly, double semiMajorAxis, double eccentricity, double inclination, double longitudeOfAscendingNode, double argumentOfPeriapsis)
        {
            this.name = name;

            this.Mass = Mass;

            this.startingTimeFromEpoch = startingTimeFromEpoch;

            failed = InitialiseOrbitConstantsFromKeplerianElements(tree, ParentMass, parentHillSphereRadius, parentPosition, startingTrueAnomaly, semiMajorAxis, eccentricity, inclination, longitudeOfAscendingNode, argumentOfPeriapsis);
        }

        private bool InitialiseOrbitConstantsFromStateVectors(Tree tree, double ParentMass, double parentHillsSphereRadius, Vector3 startingRelativePosition, Vector3 startingRelativeVelocity, Vector3 startingAbsolutePosition)
        {
            double startingRelativeAbsPosition = startingRelativePosition.Length();

            double startingRelativeAbsVelocity = startingRelativeVelocity.Length();

            // This is the vector perpendicular from the other 2 vectors which is a measure of angular momentum / mass.
            Vector3 SpecificAngularMomentum = Vector3.Cross(startingRelativePosition, startingRelativeVelocity);

            // e is the eccentricity vector, a vector with no unit that points from apoapsis to periapsis with magnitude equal to scalar eccentricity.
            Vector3 eccentricityVector = Vector3.Divide(Vector3.Cross(startingRelativeVelocity, SpecificAngularMomentum), Convert.ToSingle(G * ParentMass)) - Vector3.Divide(startingRelativePosition, Convert.ToSingle(startingRelativeAbsPosition));

            // Orbital eccentricity is a unitless (dimensionless) parameter which gives a measure of how far an orbit is from a perfect circle.
            double eccentricity = eccentricityVector.Length();

            // This is the distance from the centre of the orbit and the apoapsis or periapsis (the longest diameter).
            double semiMajorAxis = 1 / (2 / startingRelativeAbsPosition - Math.Pow(startingRelativeAbsVelocity, 2) / (G * ParentMass));

            //** Insert function for SoI calculation here;
            double hillSphereRadius = semiMajorAxis * (1 - eccentricity) * Math.Pow(Mass / (3 * ParentMass), (double)1 / 3);

            // Check that no other bodies already exist in that hill sphere.
            if (tree.AreObjectsWithinRadius(hillSphereRadius, startingAbsolutePosition))
            {
                return true;
            }

            // This is the shortest distance from the parent body.
            // **There is an acceptable error in this calculation added on to stop edge case calculation errors caused by rounding errors.
            double periapsis = semiMajorAxis * (1 - eccentricity) + 0.0001;

            // This is the largest distance from the parent body.
            // **There is an acceptable error in this calculation added on to stop edge case calculation errors caused by rounding errors.
            double apoapsis = semiMajorAxis * (1 + eccentricity) - 0.0001;

            if (apoapsis + hillSphereRadius > parentHillsSphereRadius)
            {
                //** Warn that any subsequently made bodies must lie in the parent’s hill sphere before being accepted into and child’s gravitational influence.
            }

            // The semi-minor axis is the shortest diameter of an ellipse. It is perpendicular to the semi-major axis through the centre point if drawn on a graph.
            double semiMinorAxis = semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity, 2));

            // The angle between the equatorial plane (in this case xy) and the plane of the orbit.
            double inclination = Math.Acos(SpecificAngularMomentum.Z / SpecificAngularMomentum.Length());

            // n is the vector pointing towards the ascending node. If the ascending node doesn't exist (inclination = 0 or pi) it is assigned as the nullVector3. It is used for the following equations.
            Vector3 n = new Vector3();
            if (inclination != 0 && inclination != Math.PI)
            {
                n = Vector3.Cross(new Vector3(0, 0, 1), SpecificAngularMomentum);
            }

            // The angle between the +x direction and the ascending node.
            double longitudeOfAscendingNode = 0;
            if (inclination != 0 && inclination != Math.PI)
            {
                if (n.Y >= 0)
                {
                    longitudeOfAscendingNode = Math.Acos(n.X / n.Length());
                }
                else
                {
                    longitudeOfAscendingNode = 2 * Math.PI - Math.Acos(n.X / n.Length());
                }
            }


            // The angle between the +x direction and the periapsis. This is being used as if inclination == 0 the argument of periapsis defaults to 0.
            double argumentOfPeriapsis;
            if (inclination == 0 || inclination == Math.PI)
            {
                argumentOfPeriapsis = Math.Atan2(eccentricityVector.Y, eccentricityVector.X);

                if (Vector3.Cross(startingRelativePosition, startingRelativeVelocity).Z < 0)
                {
                    argumentOfPeriapsis = 2 * Math.PI - argumentOfPeriapsis;
                }
            }
            else
            {
                argumentOfPeriapsis = Math.Acos(Vector3.Dot(n, eccentricityVector) / (n.Length() * eccentricityVector.Length()));

                if (eccentricityVector.Z < 0)
                {
                    argumentOfPeriapsis = 2 * Math.PI - argumentOfPeriapsis;
                }
            }

            double longitudeOfPeriapsis = longitudeOfAscendingNode + argumentOfPeriapsis;

            double startingPerpendicularAbsVelocity;
            if (inclination != 0)
            {
                // This method quite possibly is overcomplicating the maths. I'll need to revist this later.
                Vector3 flattenedRelativeStartPosition = startingRelativePosition;
                Vector3 flattenedRelativeStartVelocity = startingRelativeVelocity;

                // This rotates the line between the ascending and decending node to be along the positive x axis.
                double angleFromAscendingNode = 2 * Math.PI - longitudeOfAscendingNode;
                Matrix3x3 rotationMatrixForAscendingNode = new Matrix3x3(new double[,] { { Math.Cos(angleFromAscendingNode), Math.Sin(angleFromAscendingNode), 0 }, { -Math.Sin(angleFromAscendingNode), Math.Cos(angleFromAscendingNode), 0 }, { 0, 0, 1 } });

                flattenedRelativeStartPosition = rotationMatrixForAscendingNode.Transform(flattenedRelativeStartPosition);
                flattenedRelativeStartVelocity = rotationMatrixForAscendingNode.Transform(flattenedRelativeStartVelocity);

                // This rotates the orbit along it's node line (the x axis), hence giving a flat orbit.
                double angleFromInclination = 2 * Math.PI - inclination;
                Matrix3x3 rotationMatrixForInclination = new Matrix3x3(new double[,] { { 1, 0, 0 }, { 0, Math.Cos(angleFromInclination), Math.Sin(angleFromInclination) }, { 0, -Math.Sin(angleFromInclination), Math.Cos(angleFromInclination) } });

                flattenedRelativeStartPosition = rotationMatrixForInclination.Transform(flattenedRelativeStartPosition);
                flattenedRelativeStartVelocity = rotationMatrixForInclination.Transform(flattenedRelativeStartVelocity);

                // This rotates the orbit such that its periapsis and apoapsis lie on the x axis, specifically such that the periapsis is furthest along the positive x direction.
                double angleFromPeriapsis = 2 * Math.PI - argumentOfPeriapsis;
                Matrix3x3 rotationMatrixForPeriapsis = new Matrix3x3(new double[,] { { Math.Cos(angleFromPeriapsis), Math.Sin(angleFromPeriapsis), 0 }, { -Math.Sin(angleFromPeriapsis), Math.Cos(angleFromPeriapsis), 0 }, { 0, 0, 1 } });

                flattenedRelativeStartPosition = rotationMatrixForPeriapsis.Transform(flattenedRelativeStartPosition);
                flattenedRelativeStartVelocity = rotationMatrixForPeriapsis.Transform(flattenedRelativeStartVelocity);

                //startingAngleOfSpeed = Math.PI / 2 + Math.Atan(flattenedRelativeStartPosition.Y / flattenedRelativeStartPosition.X) - Math.Atan(flattenedRelativeStartVelocity.Y / flattenedRelativeStartVelocity.X);

                double a = flattenedRelativeStartPosition.Y / flattenedRelativeStartPosition.X;
                double b = flattenedRelativeStartVelocity.Y / flattenedRelativeStartVelocity.X;

                double cosStartingAngleOfSpeed = (a / Math.Sqrt(Math.Pow(a, 2) + 1)) * (1 / Math.Sqrt(Math.Pow(b, 2) + 1)) - (1 / Math.Sqrt(Math.Pow(a, 2) + 1)) * (b / Math.Sqrt(Math.Pow(b, 2) + 1));

                //startingPerpendicularAbsVelocity = startingRelativeAbsVelocity * Math.Cos(startingAngleOfSpeed);
                startingPerpendicularAbsVelocity = startingRelativeAbsVelocity * cosStartingAngleOfSpeed;
            }
            else
            {
                //startingAngleOfSpeed = Math.PI / 2 + Math.Atan(startingRelativePosition.Y / startingRelativePosition.X) - Math.Atan(startingRelativeVelocity.Y / startingRelativeVelocity.X);

                double a = startingRelativePosition.Y / startingRelativePosition.X;
                double b = startingRelativeVelocity.Y / startingRelativeVelocity.X;

                double cosStartingAngleOfVelocity = (a / Math.Sqrt(Math.Pow(a, 2) + 1)) * (1 / Math.Sqrt(Math.Pow(b, 2) + 1)) - (1 / Math.Sqrt(Math.Pow(a, 2) + 1)) * (b / Math.Sqrt(Math.Pow(b, 2) + 1));

                //startingPerpendicularAbsVelocity = startingRelativeAbsVelocity * Math.Cos(startingAngleOfSpeed);
                startingPerpendicularAbsVelocity = startingRelativeAbsVelocity * cosStartingAngleOfVelocity;
            }

            double AngularMomentum = Mass * startingPerpendicularAbsVelocity * startingRelativeAbsPosition;

            double TotalEnergy = 0.5 * Mass * Math.Pow(startingRelativeAbsVelocity, 2) + -1 * G * ParentMass * Mass / startingRelativeAbsPosition;

            currentPoint = new PointInformation(time: startingTimeFromEpoch, relativePosition: startingRelativePosition, relativeVelocity: startingRelativeVelocity);

            double orbitalPeriod;

            predeterminedPoints = createPoints(new OrbitInformation(longitudeOfAscendingNode, argumentOfPeriapsis, longitudeOfPeriapsis, inclination, eccentricity, semiMajorAxis, semiMinorAxis, periapsis, apoapsis, AngularMomentum, SpecificAngularMomentum, eccentricityVector, TotalEnergy, hillSphereRadius, 0), startingTimeFromEpoch, 20000, ParentMass, startingRelativeAbsPosition, out orbitalPeriod);

            orbitInformation = new OrbitInformation(longitudeOfAscendingNode, argumentOfPeriapsis, longitudeOfPeriapsis, inclination, eccentricity, semiMajorAxis, semiMinorAxis, periapsis, apoapsis, AngularMomentum, SpecificAngularMomentum, eccentricityVector, TotalEnergy, hillSphereRadius, (double)orbitalPeriod);

            return false;
        }

        private bool InitialiseOrbitConstantsFromKeplerianElements(Tree tree, double ParentMass, double parentHillsSphereRadius, Vector3 parentPosition, double startingTrueAnomaly, double semiMajorAxis, double eccentricity, double inclination, double longitudeOfAscendingNode, double argumentOfPeriapsis)
        {
            double startingRelativeAbsPosition = (semiMajorAxis * (1 - Math.Pow(eccentricity, 2))) / (1 + eccentricity * Math.Cos(startingTrueAnomaly));
            Vector3 startingRelativePosition = new Vector3(Convert.ToSingle(startingRelativeAbsPosition * Math.Cos(startingTrueAnomaly)), Convert.ToSingle(startingRelativeAbsPosition * Math.Sin(startingTrueAnomaly)), 0);

            // This is the shortest distance from the parent body.
            // **There is an acceptable error in this calculation added on to stop edge case calculation errors caused by rounding errors.
            double apoapsis = semiMajorAxis * (1 + eccentricity);
            //apoapsis -= apoapsis / 40000;
            double periapsis = semiMajorAxis * (1 - eccentricity);
            //periapsis += periapsis / 40000;

            Vector3 periapsisRelativePosition = new Vector3(Convert.ToSingle(periapsis), 0, 0);

            double totalVelocitySquaredAtPeriapsis = G * ParentMass * (2 / periapsis - 1 / semiMajorAxis);

            double TotalEnergy = (-G * ParentMass * Mass) / (apoapsis + periapsis);

            Vector3 periapsisRelativeVelocity = new Vector3(0, Convert.ToSingle(Math.Sqrt(totalVelocitySquaredAtPeriapsis)), 0);

            Matrix3x3 rotationMatrixForPeriapsis = new Matrix3x3(new double[,] { { Math.Cos(argumentOfPeriapsis), Math.Sin(argumentOfPeriapsis), 0 }, { -Math.Sin(argumentOfPeriapsis), Math.Cos(argumentOfPeriapsis), 0 }, { 0, 0, 1 } });

            periapsisRelativePosition = rotationMatrixForPeriapsis.Transform(periapsisRelativePosition);
            periapsisRelativeVelocity = rotationMatrixForPeriapsis.Transform(periapsisRelativeVelocity);

            Matrix3x3 rotationMatrixForInclination = new Matrix3x3(new double[,] { { 1, 0, 0 }, { 0, Math.Cos(inclination), Math.Sin(inclination) }, { 0, -Math.Sin(inclination), Math.Cos(inclination) } });

            periapsisRelativePosition = rotationMatrixForInclination.Transform(periapsisRelativePosition);
            periapsisRelativeVelocity = rotationMatrixForInclination.Transform(periapsisRelativeVelocity);

            Matrix3x3 rotationMatrixForAscendingNode = new Matrix3x3(new double[,] { { Math.Cos(longitudeOfAscendingNode), Math.Sin(longitudeOfAscendingNode), 0 }, { -Math.Sin(longitudeOfAscendingNode), Math.Cos(longitudeOfAscendingNode), 0 }, { 0, 0, 1 } });

            periapsisRelativePosition = rotationMatrixForAscendingNode.Transform(periapsisRelativePosition);
            periapsisRelativeVelocity = rotationMatrixForAscendingNode.Transform(periapsisRelativeVelocity);

            //** Insert function for SoI calculation here;
            double hillSphereRadius = semiMajorAxis * (1 - eccentricity) * Math.Pow(Mass / (3 * ParentMass), (double)1 / 3);

            Vector3 startingAbsolutePosition = startingRelativePosition + parentPosition;

            // Check that no other bodies already exist in that hill sphere.
            if (tree.AreObjectsWithinRadius(hillSphereRadius, startingAbsolutePosition))
            {
                return true;
            }

            if (apoapsis + hillSphereRadius > parentHillsSphereRadius)
            {
                //** Warn that any subsequently made bodies must lie in the parent’s hill sphere before being accepted into and child’s gravitational influence.
            }

            // The semi-minor axis is the shortest diameter of an ellipse. It is perpendicular to the semi-major axis through the centre point if drawn on a graph.
            double semiMinorAxis = semiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity, 2));

            // This is the vector perpendicular from the other 2 vectors which is a measure of angular momentum / mass.
            Vector3 SpecificAngularMomentum = Vector3.Cross(periapsisRelativePosition, periapsisRelativeVelocity);

            // e is the eccentricity vector, a vector with no unit that points from apoapsis to periapsis with magnitude equal to scalar eccentricity.
            Vector3 eccentricityVector = Vector3.Divide(Vector3.Cross(periapsisRelativeVelocity, SpecificAngularMomentum), Convert.ToSingle(G * ParentMass)) - Vector3.Divide(periapsisRelativePosition, periapsisRelativePosition.Length());

            double longitudeOfPeriapsis = longitudeOfAscendingNode + argumentOfPeriapsis;

            double startingPeriapsisAbsVelocity = periapsisRelativeVelocity.Length();

            double AngularMomentum = Mass * startingPeriapsisAbsVelocity * periapsisRelativePosition.Length();

            double startingPerpendicularVelocity = AngularMomentum / (Mass * startingRelativeAbsPosition);
            double startingRelativeAbsVelocitySquared = G * ParentMass * (2 / periapsis - 1 / semiMajorAxis);

            if (startingRelativeAbsVelocitySquared <= periapsis)
            {
                startingPerpendicularVelocity += startingPerpendicularVelocity / 400000;
            }
            else if (startingRelativeAbsVelocitySquared >= apoapsis)
            {
                startingPerpendicularVelocity -= startingPerpendicularVelocity / 400000;
            }

            Vector3 startingRelativeVelocity = new Vector3(Convert.ToSingle(startingPerpendicularVelocity), Convert.ToSingle(Math.Sqrt(Math.Pow(startingPerpendicularVelocity, 2) - startingRelativeAbsVelocitySquared)), 0);

            if (startingTrueAnomaly > Math.PI || startingTrueAnomaly < 0)
            {
                startingRelativeVelocity *= -1;
            }

            startingRelativeVelocity = rotationMatrixForPeriapsis.Transform(startingRelativeVelocity);
            startingRelativeVelocity = rotationMatrixForInclination.Transform(startingRelativeVelocity);
            startingRelativeVelocity = rotationMatrixForAscendingNode.Transform(startingRelativeVelocity);

            currentPoint = new PointInformation(time: startingTimeFromEpoch, relativePosition: startingRelativePosition, relativeVelocity: startingRelativeVelocity);

            double orbitalPeriod;

            predeterminedPoints = createPoints(new OrbitInformation(longitudeOfAscendingNode, argumentOfPeriapsis, longitudeOfPeriapsis, inclination, eccentricity, semiMajorAxis, semiMinorAxis, periapsis, apoapsis, AngularMomentum, SpecificAngularMomentum , eccentricityVector, TotalEnergy, hillSphereRadius, 0), startingTimeFromEpoch, 200000, ParentMass, startingRelativeAbsPosition, out orbitalPeriod);

            orbitInformation = new OrbitInformation(longitudeOfAscendingNode, argumentOfPeriapsis, longitudeOfPeriapsis, inclination, eccentricity, semiMajorAxis, semiMinorAxis, periapsis, apoapsis, AngularMomentum, SpecificAngularMomentum, eccentricityVector, TotalEnergy, hillSphereRadius, (double)orbitalPeriod);

            return false;
        }

        private PointInformation[] createPoints(OrbitInformation orbitInformation, int timeOffset, int numberOfIntervals, double ParentMass, double startingDistance, out double orbitalPeriod)
        {
            // V is the specific energy.
            double V = orbitInformation.TotalEnergy / Mass;

            // A is the absolute value of the specific angular momentum.
            // It is derived using -2 * V * U / K = 0 (bare in mind U is the only important part as it contains A).
            // This is because of the "time" equation used later, which requires virtually 0 error in A, breaks when using the derivation from A = Angular Momentum / Mass due to the error being increased.
            // The area that breaks is Math.Sqrt(-2 * V * U / K)  hence needs to equal basically 0 as to stretch the time equation's bounds to be as close as possible to correct.
            double A = Math.Sqrt(2 * V * Math.Pow(orbitInformation.periapsis, 2) + 2 * G * ParentMass * orbitInformation.periapsis);

            if (A > 0)
            {
                A -= A * Math.Pow(10, -16);
            }
            else
            {
                A += A * Math.Pow(10, -16);
            }

            // K is a constant used in the following equations.
            double K = Math.Pow(G * ParentMass, 2) + 2 * V * Math.Pow(A, 2);

            // midpointDistance is the distance halfway between apoapsis and periapsis.
            double midpointDistance = -G * ParentMass / (V * 2);

            // U is a variable determined by r; midpointAdjustedU is hence a variable determined by midpointDistance.
            double midpointAdjustedU = 2 * V * Math.Pow(midpointDistance, 2) + 2 * G * ParentMass * midpointDistance - Math.Pow(A, 2);

            // U is a variable determined by r; periapsisU is hence a variable determined by periapsis.
            double periapsisU = 2 * V * Math.Pow(orbitInformation.periapsis - orbitInformation.periapsis / 400000, 2) + 2 * G * ParentMass * orbitInformation.periapsis - Math.Pow(A, 2);

            // midpointTimeAdjustment is the time added to the second equation for time too adjust for the discontinuity in the graphs.
            double midpointTimeAdjustment = -2 * G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * midpointAdjustedU / K), 5)));

            // orbitalPeriod is the total time taken to complete one orbit.
            orbitalPeriod = 2 * Math.PI * Math.Sqrt(Math.Pow(orbitInformation.semiMajorAxis, 3) / (G * ParentMass));
            //**orbitalPeriod = 2 * (G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * periapsisU / K), 7))) + 1 / (2 * V) * Math.Sqrt(periapsisU) + midpointTimeAdjustment);

            /*decimal U1 = 2 * (decimal)V * (decimal)orbitInformation.periapsis * (decimal)orbitInformation.periapsis + 2 * (decimal)G * (decimal)ParentMass * (decimal)orbitInformation.periapsis - (decimal)Math.Pow(A, 2);
            decimal periapsisTime = (-1 * (decimal)G * (decimal)ParentMass / (4 * (decimal)V) * ((decimal)Math.Sqrt(-2 / V) * (decimal)Math.Asin(Math.Round(Math.Sqrt(-2 * V * (double)U1 / K), 5))) + 1 / (2 * (decimal)V) * (decimal)Math.Sqrt((double)U1) + timeOffset);

            decimal U2 = 2 * (decimal)V * (decimal)orbitInformation.apoapsis * (decimal)orbitInformation.apoapsis + 2 * (decimal)G * (decimal)ParentMass * (decimal)orbitInformation.apoapsis - (decimal)Math.Pow(A, 2);
            decimal variable = ((decimal)Math.Asin(Math.Sqrt(-2 * V * (double)U2 / K)));
            decimal apoapsisTime = ((decimal)G * (decimal)ParentMass / (4 * (decimal)V) * ((decimal)Math.Sqrt(-2 / V)) * variable + 1 / (2 * (decimal)V) * (decimal)Math.Sqrt((double)U2) + (decimal)midpointTimeAdjustment + timeOffset);

            orbitalPeriod = (apoapsisTime - periapsisTime) * 2;*/

            double range = orbitInformation.apoapsis - orbitInformation.periapsis;

            // Given that the position the body starts at is not the apoapsis, there is an offset in the time from the natural time given from starting at apoapsis.
            double startPositionTimeOffset = 0;
            // If at the extremes errors occur right on those edges. This means a small offset is given to keep the edgecases from throwing errors.
            if (startingDistance <= orbitInformation.periapsis)
            {
                CalculatePoint(startingDistance + startingDistance / 400000, startPositionTimeOffset, orbitalPeriod, true, out startPositionTimeOffset);
            }
            else if (startingDistance >= orbitInformation.apoapsis)
            {
                CalculatePoint(startingDistance - startingDistance / 400000, startPositionTimeOffset, orbitalPeriod, true, out startPositionTimeOffset);
            }
            else
            {
                CalculatePoint(startingDistance, startPositionTimeOffset, orbitalPeriod, true, out startPositionTimeOffset);
            }

            // When at a non-inclined orbit, standard methods of calculating the argument of periapsis do not work due to their reliance on the ascending node's position.
            // A different method can be used; given a starting position you can calculate it's angle from the positive x direction as well as its true anomaly, hence finding the difference to calculate the overall offset.
            // Note that strictly speaking given the definition of the AoP, this angle is to the ascending node; in this situation it is a psuedo-node for there is no inclination.

            PointInformation[] predeterminedPoints = new PointInformation[numberOfIntervals];
            /*for (int i = 0; i < numberOfIntervals; i += 2)
            {
                // Points lie both above and below the semi-major axis line, hence a positive and negative case need to be created.
                double R = orbitInformation.periapsis + i * range / numberOfIntervals;

                double positiveTime;
                double negativeTime;

                CalculatePoint(R, startPositionTimeOffset, orbitalPeriod, true, out positiveTime);
                CalculatePoint(R, startPositionTimeOffset, orbitalPeriod, false, out negativeTime);

                /*if (Math.Abs(positiveTime) < 5 || Math.Abs(negativeTime) < 5)
                {
                    // Inclination is flipping the start location because it is.
                }

                double angleOfVelocityToTangent = Math.Acos(orbitInformation.AngularMomentum / Math.Sqrt(2 * orbitInformation.TotalEnergy * Mass * Math.Pow(R, 2) + 2 * G * ParentMass * Math.Pow(Mass, 2) * R));

                double positiveTrueAnomaly = Math.Acos((orbitInformation.semiMajorAxis * (1 - Math.Pow(orbitInformation.eccentricity, 2)) - R) / (orbitInformation.eccentricity * R));
                double negativeTrueAnomaly = 2 * Math.PI - positiveTrueAnomaly;

                Vector3 positiveRelativePosition = new Vector3(Convert.ToSingle(R * Math.Cos(positiveTrueAnomaly)), Convert.ToSingle(R * Math.Sin(positiveTrueAnomaly)), 0);
                Vector3 negativeRelativePosition = new Vector3(Convert.ToSingle(R * Math.Cos(negativeTrueAnomaly)), Convert.ToSingle(R * Math.Sin(negativeTrueAnomaly)), 0);

                double perpendicularVelocity = orbitInformation.AngularMomentum / (Mass * R);
                double totalVelocitySquared = 2 * orbitInformation.TotalEnergy / Mass + (2 * G * ParentMass) / R;
                Vector3 positiveRelativeVelocity = new Vector3(Convert.ToSingle(perpendicularVelocity), Convert.ToSingle(Math.Sqrt(totalVelocitySquared - Math.Pow(perpendicularVelocity, 2))), 0);
                Vector3 negativeRelativeVelocity = new Vector3(Convert.ToSingle(-perpendicularVelocity), -Convert.ToSingle(Math.Sqrt(totalVelocitySquared - Math.Pow(perpendicularVelocity, 2))), 0);

                RotatePoint(ref positiveRelativePosition, ref positiveRelativeVelocity);
                RotatePoint(ref negativeRelativePosition, ref negativeRelativeVelocity);

                predeterminedPoints[i] = new PointInformation(positiveTime, angleOfVelocityToTangent, positiveTrueAnomaly, positiveRelativePosition, positiveRelativeVelocity);
                predeterminedPoints[i + 1] = new PointInformation(negativeTime, angleOfVelocityToTangent, negativeTrueAnomaly, negativeRelativePosition, negativeRelativeVelocity);
            }*/

            //** string output = "";
            for (int i = 0; i < numberOfIntervals; i++)
            {
                // Points lie both above and below the semi-major axis line, hence a positive and negative case need to be created.
                double trueAnomaly = i * 2 * Math.PI / numberOfIntervals;

                double relativeAbsPosition = (orbitInformation.semiMajorAxis * (1 - Math.Pow(orbitInformation.eccentricity, 2))) / (1 + orbitInformation.eccentricity * Math.Cos(trueAnomaly));
                Vector3 relativePosition = new Vector3(Convert.ToSingle(relativeAbsPosition * Math.Cos(trueAnomaly)), Convert.ToSingle(relativeAbsPosition * Math.Sin(trueAnomaly)), 0);

                double time;

                bool positive = true;

                if (trueAnomaly > Math.PI)
                {
                    positive = false;
                }

                CalculatePoint(relativeAbsPosition, startPositionTimeOffset, orbitalPeriod, positive, out time);

                double angleOfVelocityToTangent = Math.Acos(orbitInformation.AngularMomentum / Math.Sqrt(2 * orbitInformation.TotalEnergy * Mass * Math.Pow(relativeAbsPosition, 2) + 2 * G * ParentMass * Math.Pow(Mass, 2) * relativeAbsPosition));

                double perpendicularVelocity = orbitInformation.AngularMomentum / (Mass * relativeAbsPosition);
                double totalVelocitySquared = 2 * orbitInformation.TotalEnergy / Mass + (2 * G * ParentMass) / relativeAbsPosition;
                Vector3 relativeVelocity = new Vector3(Convert.ToSingle(perpendicularVelocity), Convert.ToSingle(Math.Sqrt(totalVelocitySquared - Math.Pow(perpendicularVelocity, 2))), 0);

                RotatePoint(ref relativePosition, ref relativeVelocity);

                //**
                //output += trueAnomaly + ",";
                //output += relativeAbsPosition + ",";
                //output += relativePosition + ",";
                //output += time + ",";
                //output += angleOfVelocityToTangent + ",";
                //output += relativeVelocity + "\n";

                predeterminedPoints[i] = new PointInformation(time, angleOfVelocityToTangent, trueAnomaly, relativePosition, relativeVelocity);
            }

            //** File.WriteAllText("output.txt", output);

            // **This is unnecessary, good for debugging though.
            predeterminedPoints = predeterminedPoints.OrderBy(new Func<PointInformation, double>(x => x.time)).ToArray();


            return predeterminedPoints;

            // Repeating times at 9998 and 10002?
            void CalculatePoint(double R, double startPositionTimeOffset, double orbitalPeriod, bool positive, out double time)
            {
                double U = 2 * V * R * R + 2 * G * ParentMass * R - A * A;

                //-2 * V * (2 * V * R * R + 2 * G * ParentMass * R - A * A) / K

                /*if (R < midpointDistance)
                {
                    time = (-1 * G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + timeOffset + orbitalPeriod + Math.Pow(-1, Convert.ToInt16(!positive)) * startPositionTimeOffset) % orbitalPeriod;
                }
                else
                {
                    time = (G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + midpointTimeAdjustment + timeOffset + orbitalPeriod + Math.Pow(-1, Convert.ToInt16(!positive)) * startPositionTimeOffset) % orbitalPeriod;
                }*/

                if (R < midpointDistance && positive)
                {
                    time = (-G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + timeOffset + orbitalPeriod - startPositionTimeOffset) % orbitalPeriod;
                }
                else if (R > midpointDistance && positive)
                {
                    time = (G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + midpointTimeAdjustment + timeOffset + orbitalPeriod - startPositionTimeOffset) % orbitalPeriod;
                }
                else if (R < midpointDistance && !positive)
                {
                    time = orbitalPeriod - (-G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + timeOffset + orbitalPeriod + startPositionTimeOffset) % orbitalPeriod;
                }
                else if (R > midpointDistance && !positive)
                {
                    time = orbitalPeriod - (G * ParentMass / (4 * V) * (Math.Sqrt(-2 / V) * Math.Asin(Math.Round(Math.Sqrt(-2 * V * U / K), 5))) + 1 / (2 * V) * Math.Sqrt(U) + midpointTimeAdjustment + timeOffset + orbitalPeriod + startPositionTimeOffset) % orbitalPeriod;
                }
                else
                {
                    time = double.NaN;
                }
            }

            void RotatePoint(ref Vector3 relativePosition, ref Vector3 relativeVelocity)
            {
                // This rotates the orbit such that its periapsis is at the correct angle to the positive x direction.
                Matrix3x3 rotationMatrixForPeriapsis = new Matrix3x3(new double[,] { { Math.Cos(orbitInformation.argumentOfPeriapsis), Math.Sin(orbitInformation.argumentOfPeriapsis), 0 }, { -Math.Sin(orbitInformation.argumentOfPeriapsis), Math.Cos(orbitInformation.argumentOfPeriapsis), 0 }, { 0, 0, 1 } });

                relativePosition = rotationMatrixForPeriapsis.Transform(relativePosition);
                relativeVelocity = rotationMatrixForPeriapsis.Transform(relativeVelocity);

                // This rotates the orbit to be the correct inclination.
                Matrix3x3 rotationMatrixForInclination = new Matrix3x3(new double[,] { { 1, 0, 0 }, { 0, Math.Cos(orbitInformation.inclination), Math.Sin(orbitInformation.inclination) }, { 0, -Math.Sin(orbitInformation.inclination), Math.Cos(orbitInformation.inclination) } });

                relativePosition = rotationMatrixForInclination.Transform(relativePosition);
                relativeVelocity = rotationMatrixForInclination.Transform(relativeVelocity);

                // This rotates the orbit such that the ascending node is at the correct angle to the positive x direction.
                Matrix3x3 rotationMatrixForAscendingNode = new Matrix3x3(new double[,] { { Math.Cos(orbitInformation.longitudeOfAscendingNode), Math.Sin(orbitInformation.longitudeOfAscendingNode), 0 }, { -Math.Sin(orbitInformation.longitudeOfAscendingNode), Math.Cos(orbitInformation.longitudeOfAscendingNode), 0 }, { 0, 0, 1 } });

                relativePosition = rotationMatrixForAscendingNode.Transform(relativePosition);
                relativeVelocity = rotationMatrixForAscendingNode.Transform(relativeVelocity);
            }
        }

        public void updateVariables(double time)
        {
            double[] predeterminedAdjustedTimes = new double[predeterminedPoints.Count()];

            time = (time - startingTimeFromEpoch) % orbitInformation.orbitalPeriod;

            int count = 0;
            foreach (PointInformation point in predeterminedPoints)
            {
                predeterminedAdjustedTimes[count] = Math.Abs(Convert.ToInt64(point.time) - time);
                count++;
            }

            currentPoint = predeterminedPoints[Array.IndexOf(predeterminedAdjustedTimes, predeterminedAdjustedTimes.Min())];
        }
    }
}
