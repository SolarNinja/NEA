using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Security.AccessControl;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleTestingArea
{
    internal struct Tree
    {
        private struct Node
        {
            public List<Node> childNodes { readonly get; private set; }
            //public IList<Node> childNodes { get { return new ReadOnlyCollection<Node>(_childNodes); } set { childNodes = value; } }
            public Body assignedBody { get; private set; }

            public Node(Body body)
            {
                childNodes = new List<Node>();
                assignedBody = body;
            }

            public void AddChild(Node newNode)
            {
                childNodes.Add(newNode);
            }
        }

        private Node ReferenceBody;
        public Tree(FixedBody referenceBody)
        {
            ReferenceBody = new Node(referenceBody);
        }

        public bool AreObjectsWithinRadius(double radius, Vector3 currentAbsolutePosition)
        {
            Queue<KeyValuePair<Node,Vector3>> queue = new Queue<KeyValuePair<Node,Vector3>>();
            List<Node> visitedNodes = new List<Node>();

            queue.Enqueue(new KeyValuePair<Node,Vector3>(ReferenceBody, ((FixedBody)ReferenceBody.assignedBody).position));
            visitedNodes.Add(ReferenceBody);
            while (queue.Count > 0)
            {
                KeyValuePair<Node, Vector3> nodeAndPosition = queue.Dequeue();
                //** remove this: Console.WriteLine(nodeAndPosition.Key.assignedBody.name + nodeAndPosition.Value);

                foreach (Node child in nodeAndPosition.Key.childNodes)
                {
                    if (!visitedNodes.Contains(child))
                    {
                        Vector3 childAbsolutePosition = ((MovingBody)child.assignedBody).currentPoint.relativePosition + nodeAndPosition.Value;
                        Vector3 RelativeVector = childAbsolutePosition - currentAbsolutePosition;
                        if (RelativeVector.Length() < radius)
                        {
                            return true;
                        }

                        queue.Enqueue(new KeyValuePair<Node, Vector3>(child, childAbsolutePosition));
                        visitedNodes.Add(child);
                    }
                }
            }
            return false;
        }

        public void AddToTree(string newBodyName, double newBodyMass, int time, Vector3 currentRelativeVelocity, Vector3 currentAbsolutePosition)
        {
            Tree thisTree = this;
            FindAndPopulatePosition(ReferenceBody, currentAbsolutePosition);

            void FindAndPopulatePosition(Node parentNode, Vector3 currentCalculatedPosition)
            {
                List<Node> validParents = parentNode.childNodes.FindAll(x => (currentCalculatedPosition - ((MovingBody)x.assignedBody).currentPoint.relativePosition).Length() < ((MovingBody)x.assignedBody).orbitInformation.hillSphereRadius);

                if (validParents.Count > 0)
                {
                    Node validParent = validParents[0];
                    
                    if (validParents.Count > 1)
                    {
                        validParent = validParents.MaxBy(new Func<Node, double>(x => x.assignedBody.Mass));
                    }

                    currentCalculatedPosition -= ((MovingBody)validParent.assignedBody).currentPoint.relativePosition;

                    // go again on the next scope
                    FindAndPopulatePosition(validParent, currentCalculatedPosition);
                }
                else
                {
                    // if name not taken, create new node with new body.
                    if (parentNode.childNodes.FindIndex(child => child.assignedBody.name == newBodyName) == -1)
                    {
                        bool failed = false;
                        MovingBody newBody = new MovingBody(ref failed, thisTree, newBodyName, newBodyMass, parentNode.assignedBody.Mass, ((MovingBody)parentNode.assignedBody).orbitInformation.hillSphereRadius, time, currentCalculatedPosition, currentRelativeVelocity, currentAbsolutePosition);

                        if (!failed)
                        {
                            parentNode.AddChild(new Node(newBody));
                        }
                    }
                    else
                    {
                        throw new Exception();
                        // **Output error
                    }
                }
            }
        }

        public void AddToTree(List<string> heritage, Body newBody)
        {
            // Starts at a given node, and follows a path to the newBody's parent node. Once there, the parent node is given the new body as a child node.

            // Could have used node heritage instead of searching from a tailed path; chose not too as it would require passing a recursion count as to check a specific list index, or to linearly search said list for the node of same name.
            void FindAndPopulatePosition(Node parentNode, List<string> path, Body newBody)
            {
                if (path.Count > 1)
                {
                    int childNodePosition = parentNode.childNodes.FindIndex(child => child.assignedBody.name == path[1]);
                    if (childNodePosition != -1)
                    {
                        FindAndPopulatePosition(parentNode.childNodes[childNodePosition], path.Skip(1).ToList(), newBody);
                    }
                    else
                    {
                        throw new Exception();
                        // **Output error
                    }
                }
                else if (path.Count == 1)
                {
                    if (parentNode.childNodes.FindIndex(child => child.assignedBody.name == newBody.name) == -1)
                    {
                        parentNode.AddChild(new Node(newBody));
                    }
                    else
                    {
                        throw new Exception();
                        // **Output error
                    }
                }
            }

            FindAndPopulatePosition(ReferenceBody, heritage, newBody);
        }
    }
}
