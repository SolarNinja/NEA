using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace _3D_Orbital_Motion_Simulation
{
    internal class FixedBody : Body
    {
        public override string name { get; set; }

        public Vector3 position { get; private set; }

        public override double Mass { get; set; }

        public FixedBody(string name, Vector3 position, double Mass)
        {
            this.name = name;
            this.position = position;
            this.Mass = Mass;
        }
    }
}
