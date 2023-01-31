using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.AccessControl;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

namespace _3D_Orbital_Motion_Simulation
{
    internal struct Tree
    {
        private struct Node
        {
            public List<Node> childNodes;
            public Body assignedBody;
            private List<string> heritage;

            public Node(List<string> heritage, Body body)
            {
                childNodes = new List<Node>();
                assignedBody = body;
                this.heritage = heritage;
            }



            public void AddToNode()
            {

            }
        }

        private Node ReferenceBody;
        public Tree(Body referenceBody)
        {
            ReferenceBody = new Node(null, referenceBody);
        }

        public void AddToTree(List<string> heritage, Body newBody)
        {
            // Starts at a given node, and follows a path to the newBody's parent node. Once there, the parent node is given the new body as a child node.

            // Could have used node heritage instead of searching from a tailed path; chose not too as it would require passing a recursion count as to check a specific list index, or to linearly search said list for the node of same name.
            void FindAndPopulatePosition(Node parentNode, List<string> path, Body newBody)
            {
                int childNodePosition = parentNode.childNodes.FindIndex(child => child.assignedBody.Name == path[0]);
                if (childNodePosition != -1)
                {
                    FindAndPopulatePosition(parentNode.childNodes[childNodePosition], path.Skip(1).ToList(), newBody);
                }
                else if (path.Count == 0)
                {
                    if (parentNode.childNodes.FindIndex(child => child.assignedBody.Name == newBody.Name) != -1)
                    {
                        parentNode.childNodes.Add(new Node(heritage, newBody));
                    }
                    else
                    {
                        // **Output error
                    }
                }
            }

            FindAndPopulatePosition(ReferenceBody, heritage, newBody);
        }
    }
}
