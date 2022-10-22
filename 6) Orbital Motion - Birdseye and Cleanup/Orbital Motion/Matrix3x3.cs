using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;
using System.Windows.Forms;

namespace Orbital_Motion
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

            transformedPoint.X = Convert.ToSingle(point.X * (contents[0, 0] + contents[1, 0] + contents[2, 0]));
            transformedPoint.Y = Convert.ToSingle(point.Y * (contents[0, 1] + contents[1, 1] + contents[2, 1]));
            transformedPoint.Z = Convert.ToSingle(point.Z * (contents[0, 2] + contents[1, 2] + contents[2, 2]));

            return transformedPoint;
        }
    }
}
