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
            tree = new Tree(new FixedBody("Earth", new Vector3(Convert.ToSingle(2.5 * Math.Pow(10, 7)), Convert.ToSingle(2.5 * Math.Pow(10, 7)), 0), 5.972 * Math.Pow(10, 24)));

            time = 0;

            tree.AddToTree("Luna", 7.347*Math.Pow(10, 22), time, new Vector3(1500, 1501, 0), new Vector3(25000000 - 10000000, 25000000 + 10000000, 0));
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            time += 5;

            label1.Text = time.ToString();

            tree.UpdateBodies(time);
            drawParticles();
        }

        private void drawParticles()
        {
            List<KeyValuePair<Body, Vector3>> bodies = tree.GetBodiesAndPositions();

            Graphics gr = CreateGraphics();
            Pen pen = new Pen(Color.Black);

            foreach (KeyValuePair<Body, Vector3> body in bodies)
            {
                gr.DrawRectangle(pen, Convert.ToInt64((body.Value.X) / 100000), Convert.ToInt64(Height - (body.Value.Y) / 100000), 1, 1);
            }
        }

        private void Form1_MouseClick(object sender, MouseEventArgs e)
        {
            tree.debug();
        }
    }
}