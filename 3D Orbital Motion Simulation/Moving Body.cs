using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _3D_Orbital_Motion_Simulation
{
    internal class Moving_Body : Body
    {
        public override string Name { get; set; }

        Moving_Body(string name)
        {
            this.Name = name;
        }
    }
}
