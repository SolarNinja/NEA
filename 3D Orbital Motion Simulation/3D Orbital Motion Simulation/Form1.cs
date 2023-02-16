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

        // Time is an unsigned long as it has the largest value range and it starts from 0 due to time not being negative.
        ulong time;

        private void Form1_Load(object sender, EventArgs e)
        {
            double AU = 1.49597870700 * Math.Pow(10, 11);

            tree = new Tree(new FixedBody("Sol", new Vector3(5000000000000, 5000000000000, 0), 1.989 * Math.Pow(10, 30)));

            time = 0;

            tree.AddToTree(new List<string> { "Sol" }, "Mercury", 3.285 * Math.Pow(10, 23), 0, 0, 0.3871 * AU, 0.20564, 0.122278, 0.843692, 0.50824);

            tree.AddToTree(new List<string> { "Sol" }, "Venus", 4.867 * Math.Pow(10, 24), 0, 0, 0.7233 * AU, 0.00676, 0.059306, 1.338144, 0.961676);

            tree.AddToTree(new List<string> { "Sol" }, "Earth", 5.972 * Math.Pow(10, 24), 0, 0, AU, 0.0167086, 0, 3.0525809, 5.0282936);

            tree.AddToTree(new List<string> { "Sol", "Earth" }, "Luna", 7.347 * Math.Pow(10, 22), 0, 0, 0.3844 * Math.Pow(10, 9), 0.0549, 0.08979719, 0, 0);

            tree.AddToTree(new List<string> { "Sol" }, "Mars", 6.39 * Math.Pow(10, 23), 0, 0, 1.5237 * AU, 0.09337, 0.032323, 0.867603, 4.998099);

            tree.AddToTree(new List<string> { "Sol" }, "Jupiter", 1.898 * Math.Pow(10, 27), 0, 0, 5.2025 * AU, 0.04854, 0.022672, 1.750391, -1.50133);

            tree.AddToTree(new List<string> { "Sol" }, "Saturn", 5.683 * Math.Pow(10, 26), 0, 0, 9.5415 * AU, 0.05551, 0.043529, 1.983392, -0.36268);

            tree.AddToTree(new List<string> { "Sol" }, "Uranus", 8.681 * Math.Pow(10, 25), 0, 0, 19.188 * AU, 0.04686, 0.013491, 1.290846, 1.718626);

            tree.AddToTree(new List<string> { "Sol" }, "Neptune", 1.024 * Math.Pow(10, 26), 0, 0, 30.07 * AU, 0.00895, 0.030892, 2.300169, -1.48545);
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            time += 1000000;

            YearsLabel.Text = (time / (365 * 24 * 60 * 60)).ToString();
            DaysLabel.Text = (time % (365 * 24 * 60 * 60) / (24 * 60 * 60)).ToString();
            HoursLabel.Text = (time % (365 * 24 * 60 * 60) % (24 * 60 * 60) / (60 * 60)).ToString();
            MinutesLabel.Text = (time % (365 * 24 * 60 * 60) % (24 * 60 * 60) % (60 * 60) / 60).ToString();
            SecondsLabel.Text = (time % (365 * 24 * 60 * 60) % (24 * 60 * 60) % (60 * 60) % 60).ToString();


            tree.UpdateBodies(time);
            drawParticles();
        }

        private void drawParticles()
        {
            List<KeyValuePair<Body, Vector3>> bodies = tree.GetBodiesAndPositions();

            Graphics gr = CreateGraphics();
            Pen pen;
            
            Pen[] pens = new Pen[10];

            pens[0] = new Pen(Color.Gold);
            pens[1] = new Pen(Color.DarkGray);
            pens[2] = new Pen(Color.Salmon);
            pens[3] = new Pen(Color.Blue);
            pens[4] = new Pen(Color.OrangeRed);
            pens[5] = new Pen(Color.SandyBrown);
            pens[6] = new Pen(Color.Beige);
            pens[7] = new Pen(Color.Azure);
            pens[8] = new Pen(Color.Navy);
            pens[9] = new Pen(Color.Gray);

            int count = 0;
            foreach (KeyValuePair<Body, Vector3> body in bodies)
            {
                pen = pens[count];
                count++;
                gr.DrawRectangle(pen, Convert.ToInt64((body.Value.X) / 10000000000), Convert.ToInt64(Height - (body.Value.Y) / 10000000000), 1, 1);
            }
        }

        private void Form1_MouseClick(object sender, MouseEventArgs e)
        {
            tree.debug();
        }
    }
}