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

        // Metres
        static Vector2 position = new Vector2(Convert.ToSingle(2.5 * Math.Pow(10, 7) - 6779000), Convert.ToSingle(2.5 * Math.Pow(10, 7)));
        static Vector2 velocity = new Vector2(0, 7660);

        //static double[] sun = { 2.5 * Math.Pow(10,7), 2.5 * Math.Pow(10,7), 2 * Math.Pow(6, 24) };
        static double[] sun = { 2.5 * Math.Pow(10, 7), 2.5 * Math.Pow(10, 7), 5.972 * Math.Pow(10, 24) };
        // x, y, mass

        double time = 0;

        static double G = 6.6743 * Math.Pow(10, -11);
        Body body = new Body(position, velocity, 420000, sun);

        private void Form1_Click(object sender, EventArgs e)
        {

        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            Vector2 radialPosition;
            Vector2 relativeVelocity;
            double angularVelocity;

            body.returnValues(time, sun, out radialPosition, out relativeVelocity, out angularVelocity);

            time += 5;
            
            drawParticles();
        }

        private void colourAllowed(double totalEnergy)
        {
            Graphics gr = CreateGraphics();
            Pen pen = new Pen(Brushes.Red);
            for (int x = 150; x < 350; x++)
            {
                for (int y = 200; y < 300; y++)
                {
                    double gravitationalPotentialEnergy = -1 * G * sun[2] * body.getMass() / Math.Pow(Math.Pow(x * 100000 - sun[0], 2) + Math.Pow(y * 100000 - sun[1], 2), 2);

                    if (gravitationalPotentialEnergy <= totalEnergy)
                    {
                        gr.DrawRectangle(pen, x, y, 1, 1);
                    }
                }
            }
        }

        private void sketchConic()
        {
            Graphics gr = CreateGraphics();
            Pen pen = new Pen(Brushes.Orange);

            for (double theta = 0; theta <= 2 * Math.PI; theta += 0.1)
            {

            }
        }

        private void drawParticles()
        {
            Graphics gr = CreateGraphics();
            Pen pen = new Pen(Color.Black);
            gr.DrawRectangle(pen, Convert.ToInt64(body.getPosition().X) / 100000, Convert.ToInt64(body.getPosition().Y) / 100000, 1, 1);
        }
    }
}
