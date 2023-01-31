using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _3D_Orbital_Motion_Simulation
{
    internal abstract class Body
    {
        public abstract string name { get; set; }
        public abstract double Mass { get; set; }
    }
}
