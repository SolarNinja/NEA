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

        static int startingX = -5000000;
        static int startingY = 3000000;
        static int startingVx = 5000;
        static int startingVy = 7660;
        int startingR = Convert.ToInt32(Math.Sqrt((long)startingX * startingX + (long)startingY * startingY));
        
        static Vector3 relativePosition = new Vector3(startingX, startingY, 0);
        static Vector3 velocity = new Vector3(startingVx, startingVy, 0);

        static double[] sun = { 2.5 * Math.Pow(10, 7), 2.5 * Math.Pow(10, 7), 5.972 * Math.Pow(10, 24) };

        static int NullValue = Globals.NullValue;
        static Vector3 NullVector3 = Globals.NullVector3;

        double time = 0;

        //Body body = new Body(relativePosition, velocity, 420000, sun);

        Body body = new Body(sun, 420000, 0, relativePosition, velocity, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullValue, NullVector3);

        private void Form1_Click(object sender, EventArgs e)
        {
            label6.Text = Convert.ToString(body.getRelativePosition().X) + " " + Convert.ToString(body.getRelativePosition().Y);
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            time += 5;
            body.updateVariables(time);
            label1.Text = Convert.ToString(body.getRelativePosition().X / Math.Abs(body.getRelativePosition().X));
            label2.Text = Convert.ToString(body.getRelativePosition().Y / Math.Abs(body.getRelativePosition().Y));
            label3.Text = Convert.ToString(body.getRelativePosition().X);
            label4.Text = Convert.ToString(body.getRelativePosition().Y);
            label5.Text = Convert.ToString(time + " " + body.getTime());

            label7.Text = Convert.ToString(body.getRelativeVelocity().X + " " + body.getRelativeVelocity().Y);

            if (time == 5)
            {
                label6.Text = Convert.ToString(body.getRelativePosition().X) + " " + Convert.ToString(body.getRelativePosition().Y);
            }

            drawParticles();
        }


        private void drawParticles()
        {
            Graphics gr = CreateGraphics();
            Pen pen = new Pen(Color.Black);
            gr.DrawRectangle(pen, Convert.ToInt64(sun[0] / 100000),Convert.ToInt64(Height - sun[1] / 100000), 1, 1);

            gr.DrawLine(pen, Convert.ToInt32((sun[0] + startingX) / 100000), Convert.ToInt32(Height - (sun[1] + startingY) / 100000), Convert.ToInt32((sun[0] + startingX) / 100000 + startingVx), Convert.ToInt32(Height - (sun[1] + startingY) / 100000 - startingVy));
            //gr.DrawLine(pen, Convert.ToInt32((sun[0] + startingX) / 100000), Convert.ToInt32(Height - (sun[1] + startingY) / 100000), Convert.ToInt32((sun[0] + startingX) / 100000 - startingVx), Convert.ToInt32(Height - (sun[1] + startingY) / 100000 + startingVy));

            gr.DrawRectangle(pen, Convert.ToInt64((sun[0] + startingX) / 100000), Convert.ToInt64(Height - (sun[1] + startingY) / 100000), 1, 1);
            gr.DrawEllipse(pen, Convert.ToInt32((sun[0] - startingR)/ 100000), Convert.ToInt32(Height - (sun[1] - startingR) / 100000), startingR * 2 / 100000, -startingR * 2 / 100000);

            if (Vector3.Dot(body.getRelativePosition(), body.getRelativeVelocity()) < 0)
            {
                pen = new Pen(Color.Red);
            }
            else
            {
                pen = new Pen(Color.Pink);
            }
            gr.DrawRectangle(pen, Convert.ToInt64((sun[0] + body.getRelativePosition().X) / 100000), Convert.ToInt64(Height - (sun[1] + body.getRelativePosition().Y) / 100000), 1, 1);
        }
    }
}
