package toolClasses;

import java.util.Arrays;
import java.util.LinkedList;

/**
 * Adapted Neighbour Joining algorithm after:
 * N Saitou, M Nei,
 * The neighbor-joining method: a new method for reconstructing phylogenetic trees.
 * Molecular Biology and Evolution, Volume 4, Issue 4, Jul 1987, Pages 406–425
 * https://doi.org/10.1093/oxfordjournals.molbev.a040454
 */

public class NeighbourJoiningAlgorithmus {

    private final DistanceMatrix distanceMatrix;
    private final int[][] initialDistanceMatrix;
    private final Node[] initialNodes;

    public NeighbourJoiningAlgorithmus(LinkedList<Profile> initialProfiles) {
        distanceMatrix = new DistanceMatrix(initialProfiles);
        initialDistanceMatrix = distanceMatrix.distanceMatrix;
        initialNodes = distanceMatrix.getNodes();
    }

    /**
     * starts the NJ algorithm
     */
    public void runAlgorithm() {

        while (distanceMatrix.size() != 2) {
            // compute Neighbour Matrix
            NeighbourMatrix nMatrix = new NeighbourMatrix(distanceMatrix);
            // find smallest neighbour-distance between two nodes & merge them to a new Node
            Node combinedNode = nMatrix.findAndCombineNearestNodes();
            // update distanceMatrix
            distanceMatrix.updateDistanceMatrix(combinedNode);
        }

        // merge the last two remaining Nodes
        Node root = new Node(distanceMatrix.getNodes()[0], distanceMatrix.getNodes()[1]);

        System.out.println(root.name);
    }

    /**
     * for debug purposes prints the initial
     */
    public void printInitialDistanceMatrix() {
        System.out.println("Nodes with index: ");
        int index = 0;
        for (Node node : initialNodes) {
            System.out.println("#" + index++ + " " + node.name);
        }

        System.out.println("\nMatrix:");
        for (int[] row : initialDistanceMatrix) {
            System.out.println(Arrays.toString(row));
        }
    }






    private static class DistanceMatrix {
        // symmetric distance matrix
        private int[][] distanceMatrix;
        // this list is maintained to find the corresponding index of the distance list of a node
        private Node[] nodes;

        public DistanceMatrix(LinkedList<Profile> initialSequences) {
            // put each sequence into a leaf node
            nodes = new Node[initialSequences.size()];
            for (int i = 0; i < nodes.length; i++) {
                Node node = new Node(initialSequences.get(i));
                nodes[i] = node;
            }

            // instance matrix
            distanceMatrix = new int[nodes.length][nodes.length];

            computeInitialDistances(initialSequences);
        }

        /**
         * initializes the matrix with all distances
         * @param initialSequences
         */
        public void computeInitialDistances(LinkedList<Profile> initialSequences) {
            // for each Leaf, compute the distance to each other leaf
            for (int i = 0; i < nodes.length; i++) {
                for (int j = 0; j < nodes.length; j++) {
                    if (j < i) {
                        // the d matrix is symmetric, so here we just place the already computed score into
                        // the correct index
                        distanceMatrix[i][j] = distanceMatrix[j][i];

                    } else if (i == j) {
                        distanceMatrix[i][j] = 0;

                    }  else {
                        Node node1 = nodes[i];
                        Node node2 = nodes[j];
                        String sequence1 = node1.getProfile().getInitialSequence();
                        String sequence2 = node2.getProfile().getInitialSequence();
                        int score = SequenceAlignment.computeAlignmentScore(sequence1, sequence2);
                        distanceMatrix[i][j] = score;
                    }
                }
            }
        }

        /**
         * updates the current distanceMatrix with the new Node, calculates the new distances & overrides the current
         * Matrix.
         * @param newNode newNode that holds the merged nodes as childs
         */
        public void updateDistanceMatrix(Node newNode) {
            Node[] newNodes = new Node[nodes.length - 1];
            int[][] newDistanceMatrix = new int[newNodes.length][newNodes.length];

            Node mergedNode1 = newNode.getChildNode1();
            Node mergedNode2 = newNode.getChildNode2();

            int index = 0;
            for (Node oldNode : nodes) {
                if (!oldNode.equals(mergedNode1) && !oldNode.equals(mergedNode2)) {
                    newNodes[index++] = oldNode;
                }
            }

            newNodes[newNodes.length - 1] = newNode;

            for (int i = 0; i < newNodes.length; i++) {
                for (int j = i + 1; j < newNodes.length; j++) {
                    if (j < i) {
                        // symmetric filling
                        distanceMatrix[i][j] = distanceMatrix[j][i];

                    } else if (i == j) {
                        distanceMatrix[i][j] = 0;

                    }  else {
                        Node node1 = newNodes[i];
                        Node node2 = newNodes[j];

                        // if booth nodes are "old" then we just access their values
                        if (!node1.equals(newNode) && !node2.equals(newNode)) {
                            newDistanceMatrix[i][j] = distanceMatrix[getIndexOfNode(node1)][getIndexOfNode(node2)];

                        // if one of them is the new (merged) node, we update the distance
                        } else {
                            int distance = computeNewDistanceBetween(node1, node2);
                            newDistanceMatrix[i][j] = distance;
                        }
                    }
                }
            }

            // assign the new nodes and the new matrix
            nodes = newNodes;
            distanceMatrix = newDistanceMatrix;
        }

        public int size() {
            return nodes.length;
        }

        /**
         * finds the index of the given node in the node array
         * @param node which index should be found
         * @return index of node
         */
        private int getIndexOfNode(Node node) {
            int index = 0;
            for (int i = 0; i < nodes.length; i++) {
                if (nodes[i].equals(node)) index = i;
            }
            return index;
        }

        /**
         * returns the distance between two nodes
         * @param node1
         * @param node2
         * @return distance between those two nodes
         */
        public int getDistanceBetween(Node node1, Node node2) {
            int indexNode1 = 0;
            int indexNode2 = 0;

            for (int i = 0; i < nodes.length; i++) {
                if (nodes[i].equals(node1)) indexNode1 = i;
                if (nodes[i].equals(node2)) indexNode2 = i;
            }

            return distanceMatrix[indexNode1][indexNode2];
        }

        /**
         * updates the distance between an "old" node with the new merged node
         * @param node1
         * @param newNode
         * @return the distance between them
         */
        public int computeNewDistanceBetween(Node node1, Node newNode) {
            Node node2 = newNode.getChildNode1();
            Node node3 = newNode.getChildNode2();

            return (this.getDistanceBetween(node1, node2) + this.getDistanceBetween(node1, node3) - this.getDistanceBetween(node2, node3)) / 2;
        }

        /**
         * returns the distances to every node
         * @param node
         * @return
         */
        public int[] getDistancesOf(Node node) {
            int index = this.getIndexOfNode(node);
            return distanceMatrix[index];
        }

        /**
         * returns the list of current Nodes
         * @return
         */
        public Node[] getNodes() {
            return nodes;
        }
    }









    private static class NeighbourMatrix {
        // this list is maintained to find the corresponding index of the distance list of a node
        private final Node[] nodes;
        private final DistanceMatrix distanceMatrix;
        private final int[][] neighbourMatrix;

        public NeighbourMatrix(DistanceMatrix d) {
            this.nodes = d.getNodes();
            this.distanceMatrix = d;
            this.neighbourMatrix = new int[nodes.length][nodes.length];
            computeMatrix();
        }

        private void computeMatrix() {
            // for each Leaf, compute the distance to each other leaf
            for (int i = 0; i < nodes.length; i++) {
                for (int j = 0; j < nodes.length; j++) {

                    if (j < i) {
                        // the d matrix is symmetric, so here we just place the already computed score into
                        // the correct index
                        neighbourMatrix[i][j] = neighbourMatrix[j][i];

                    } else if (j == i) {
                        neighbourMatrix[i][j] = 0;

                    }  else {
                        Node node1 = nodes[i];
                        Node node2 = nodes[j];
                        neighbourMatrix[i][j] = calculateNValue(node1, node2);
                    }
                }
            }
        }

        /**
         * find the Node pair with smallest neighbour value
         * @return a new node which its childs are the two combined nodes.
         */
        public Node findAndCombineNearestNodes() {
            Node node1OfSmallestPair = null;
            Node node2OfSmallestPair = null;
            Integer smallestDistance = Integer.MAX_VALUE;
            for (int i = 0; i < nodes.length; i++) {
                for (int j = i; j < nodes.length; j++) {
                    if (neighbourMatrix[i][j] < smallestDistance) {
                        smallestDistance = neighbourMatrix[i][j];
                        node1OfSmallestPair = nodes[i];
                        node2OfSmallestPair = nodes[j];
                    }
                }
            }
            return new Node(node1OfSmallestPair, node2OfSmallestPair);
        }

        /**
         * Computes the mean distrance to every other Profile (r)
         * @param node
         * @return
         */
        private int computeDistanceMeanToEveryOtherNode(Node node) {
            int[] distances = distanceMatrix.getDistancesOf(node);
            int sum = 0;
            for (Integer distance : distances) {
                sum += distance;
            }
            return sum / (distanceMatrix.size() - 2);
        }

        /**
         * calculates the neighbour value between two nodes
         * @param node1
         * @param node2
         * @return
         */
        private int calculateNValue(Node node1, Node node2) {
            int rNode1 = computeDistanceMeanToEveryOtherNode(node1);
            int rNode2 = computeDistanceMeanToEveryOtherNode(node2);
            int distance = distanceMatrix.getDistanceBetween(node1, node2);

            return distance - (rNode1 + rNode2);
        }

    }





    // build a tree directly inside NJ (using Nodes that hold the Profiles)
    // this tree can contain a method "to Newick" for extension
    /**
     * Nodes of the guide tree
     */
    public static class Node {
        private Node parentNode;
        private Node childNode1;
        private Node childNode2;
        private Profile profile; // this is the profile this nodes holds
        private final boolean isLeaf;
        private final String name;

        public Node(Profile profile) {
            this.profile = profile;
            this.isLeaf = true;
            this.name = profile.getInitialSequence();
        }

        public Node(Node childNode1, Node childNode2) {
            this.isLeaf = false;
            this.childNode1 = childNode1;
            this.childNode2 = childNode2;
            this.name = "(" + childNode1.getName() + "," + childNode2.getName() + ")";
        }

        public String getName() {
            return name;
        }

        public boolean isLeaf() {
            return isLeaf;
        }

        public void setParentNode(Node parentNode) {
            this.parentNode = parentNode;
        }

        public void setChildNode1(Node childNode1) {
            this.childNode1 = childNode1;
        }

        public void setChildNode2(Node childNode2) {
            this.childNode2 = childNode2;
        }

        public Node getParentNode() {
            return parentNode;
        }

        public Node getChildNode1() {
            return childNode1;
        }

        public Node getChildNode2() {
            return childNode2;
        }

        public Profile getProfile() {
            return profile;
        }
    }
}
