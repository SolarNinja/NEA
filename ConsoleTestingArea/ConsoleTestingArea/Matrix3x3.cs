using System.Numerics;

namespace ConsoleTestingArea
{
    internal struct Matrix3x3
    {
        public double[,] contents;

        public Matrix3x3(double[,] contents)
        {
            if (contents.GetLength(0) != 3 || contents.GetLength(1) != 3)
            {
                throw new Exception("Incorrect dimensions inputted when creating the matrix.");
            }

            this.contents = contents;
        }

        public Vector3 Transform(Vector3 point)
        {
            Vector3 transformedPoint;

            transformedPoint.X = Convert.ToSingle(point.X * contents[0, 0] + point.Y * contents[1, 0] + point.Z * contents[2, 0]);
            transformedPoint.Y = Convert.ToSingle(point.X * contents[0, 1] + point.Y * contents[1, 1] + point.Z * contents[2, 1]);
            transformedPoint.Z = Convert.ToSingle(point.X * contents[0, 2] + point.Y * contents[1, 2] + point.Z * contents[2, 2]);

            return transformedPoint;
        }
    }
}
