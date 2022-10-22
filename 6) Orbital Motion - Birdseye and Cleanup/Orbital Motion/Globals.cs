using System;
using System.Collections.Generic;
using System.Data.SqlTypes;
using System.Numerics;
using System.Text;

namespace Orbital_Motion
{
    internal class Globals
    {
        // This value is chosen as it is mostly impossible for values to equal -1, and for those for which it can it is incredibly unlikely. This may need revision.
        static public int NullValue = -1;

        static public Vector3 NullVector3 = new Vector3(-1, -1, -1);
    }
}
