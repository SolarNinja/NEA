using System.Numerics;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.Tab;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TaskbarClock;

namespace _3D_Orbital_Motion_Simulation
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        Tree tree;
        int time;

        private void Form1_Load(object sender, EventArgs e)
        {
            tree = new Tree(new FixedBody("Sol", new Vector3(Convert.ToSingle(2.5 * Math.Pow(10, 7)), Convert.ToSingle(2.5 * Math.Pow(10, 7)), 0), 5.972 * Math.Pow(10, 24)));

            time = 0;

            tree.AddToTree("Earth", 420000, time, new Vector3(5000, 7660, 0), new Vector3(-5000000, 3000000, 0));
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            time += 5;
            tree.UpdateBodies(time);
            drawParticles();
        }

        private void drawParticles()
        {
            List<KeyValuePair<Body, Vector3>> bodies = tree.GetBodiesAndPositions();

            Graphics gr = CreateGraphics();
            Pen pen = new Pen(Color.Black);
           // gr.DrawRectangle(pen, Convert.ToInt64( / 100000), Convert.ToInt64(Height - sun[1] / 100000), 1, 1);

            //gr.DrawLine(pen, Convert.ToInt32((sun[0] + startingX) / 100000), Convert.ToInt32(Height - (sun[1] + startingY) / 100000), Convert.ToInt32((sun[0] + startingX) / 100000 + startingVx), Convert.ToInt32(Height - (sun[1] + startingY) / 100000 - startingVy));
            //gr.DrawLine(pen, Convert.ToInt32((sun[0] + startingX) / 100000), Convert.ToInt32(Height - (sun[1] + startingY) / 100000), Convert.ToInt32((sun[0] + startingX) / 100000 - startingVx), Convert.ToInt32(Height - (sun[1] + startingY) / 100000 + startingVy));

            //gr.DrawRectangle(pen, Convert.ToInt64((sun[0] + startingX) / 100000), Convert.ToInt64(Height - (sun[1] + startingY) / 100000), 1, 1);
            //gr.DrawEllipse(pen, Convert.ToInt32((sun[0] - startingR) / 100000), Convert.ToInt32(Height - (sun[1] - startingR) / 100000), startingR * 2 / 100000, -startingR * 2 / 100000);

            /*if (Vector3.Dot(body.getRelativePosition(), body.getRelativeVelocity()) < 0)
            {
                pen = new Pen(Color.Red);
            }
            else
            {
                pen = new Pen(Color.Pink);
            }*/

            foreach (KeyValuePair<Body, Vector3> body in bodies)
            {
                gr.DrawRectangle(pen, Convert.ToInt64((body.Value.X) / 100000), Convert.ToInt64(Height - (body.Value.Y) / 100000), 1, 1);
            }
        }
    }
}