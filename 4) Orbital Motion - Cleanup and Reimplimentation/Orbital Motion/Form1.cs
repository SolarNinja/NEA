using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Numerics;

namespace Orbital_Motion
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        static Vector3 relativePosition = new Vector3(6779000, 0, 0);
        static Vector3 velocity = new Vector3(0, 7660, 0);

        static double[] sun = { 2.5 * Math.Pow(10, 7), 2.5 * Math.Pow(10, 7), 5.972 * Math.Pow(10, 24) };

        static int NullValue = Globals.NullValue;
        static Vector3 NullVector3 = Globals.NullVector3;

        Body john = new Body()

        double time = 0;

        //Body body = new Body(relativePosition, velocity, 420000, sun);

        Body body = new Body(sun, relativePosition, velocity, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, 0, NullValue, NullValue, NullValue, NullValue, NullValue, 420000, NullValue, NullValue, NullValue, NullValue, NullVector3);

        private void Form1_Click(object sender, EventArgs e)
        {

        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            //Vector2 radialPosition;
            //Vector2 relativeVelocity;
            //double angularVelocity;

            //body.returnValues(time, sun);

            time += 5;
            
            drawParticles();
        }


        private void drawParticles()
        {
            Graphics gr = CreateGraphics();
            Pen pen = new Pen(Color.Black);
            //gr.DrawRectangle(pen, Convert.ToInt64(body.getPosition().X) / 100000, Convert.ToInt64(body.getPosition().Y) / 100000, 1, 1);
        }
    }
}
