package progressiveAligner.MainComponents;

import java.util.Arrays;
import java.util.LinkedList;

/**
 * Adapted Neighbour Joining algorithm after:
 * N Saitou, M Nei,
 * The neighbor-joining method: a new method for reconstructing phylogenetic trees.
 * Molecular Biology and Evolution, Volume 4, Issue 4, Jul 1987, Pages 406–425
 * <a href="https://doi.org/10.1093/oxfordjournals.molbev.a040454">...</a>
 */

public class NeighbourJoining {

    private int[][] distanceMatrix;
    private int[][] neighbourMatrix;
    private Node[] nodesOnMatrix;
    private Node root = null;
    private boolean algorithmFinished = false;

    /**
     * Constructs a new instance of the NeighbourJoining class and initializes
     * the algorithm with the given profiles.
     *
     * <p>The constructor performs the following:
     * <ul>
     *   <li>Creates an initial set of {@link Node} objects, each representing a leaf node corresponding
     *       to a profile in the input list.</li>
     *   <li>Computes the initial distance matrix by calculating pairwise distances between all profiles
     *       using sequence alignment scores.</li>
     * </ul>
     *
     * @param initialProfiles a {@link LinkedList} of {@link Profile} objects that serve as the starting
     *                        points for constructing the phylogenetic tree. Each profile corresponds
     *                        to a leaf node in the tree.
     */
    public NeighbourJoining(LinkedList<Profile> initialProfiles) {
        initialize(initialProfiles);
    }

    /**
     * Executes the Neighbour Joining (NJ) algorithm to construct a phylogenetic tree.
     *
     * <p>This method iteratively performs the following steps until only two nodes remain:
     * <ul>
     *   <li>Prints the current state of the distance matrix (for debugging purposes).</li>
     *   <li>Computes the neighbour matrix based on the current distance matrix.</li>
     *   <li>Finds the pair of nodes with the smallest neighbour-joining value and merges them into a new node.</li>
     *   <li>Updates the distance matrix to include the new merged node and recalculates distances.</li>
     * </ul>
     *
     * <p>When only two nodes remain, they are merged into a root node, completing the tree construction.
     * The resulting root node represents the entire phylogenetic tree.
     *
     * @return the root {@link Node} of the constructed phylogenetic tree.
     */
    public Node runAlgorithm() {

        // prevent second run if runAlgorithm is called twice or more.
        if (algorithmFinished) return this.root;

        while (distanceMatrix.length != 2) {
            // printCurrentDistanceMatrix(); //DEBUG
            // compute Neighbour Matrix
            computeNeighbourMatrix();
            // find smallest neighbour-distance between two nodes & merge them to a new Node
            Node combinedNode = findAndCombineNearestNodes();
            // update distanceMatrix
            updateDistanceMatrix(combinedNode);
        }

        // printCurrentDistanceMatrix(); //DEBUG

        // merge the last two remaining Nodes
        Node root = new Node(nodesOnMatrix[0], nodesOnMatrix[1]);

        // System.out.println(root.name); //DEBUG

        algorithmFinished = true;
        this.root = root;
        return root;
    }

    /**
     * for debug purposes prints the initial distance matrix
     */
    private void printCurrentDistanceMatrix() {
        System.out.println("Nodes on Matrix: ");
        int index = 0;
        for (Node node : nodesOnMatrix) {
            System.out.println("#" + index++ + " " + node.name);
        }

        System.out.println("Matrix:");
        for (int[] row : distanceMatrix) {
            System.out.println(Arrays.toString(row));
        }
        System.out.println("\n");
    }

    /**
     * Initializes the Neighbour Joining algorithm by creating nodes for the initial profiles
     * and computing the distance matrix based on the alignment scores between the sequences.
     *
     * <p>This method performs the following steps:
     * <ul>
     *   <li>Creates a {@link Node} for each {@link Profile} in the provided list of initial profiles.</li>
     *   <li>Initializes the {@code distanceMatrix} to hold pairwise distances between nodes.</li>
     *   <li>Calculates the pairwise distances based on alignment scores and populates the distance matrix.</li>
     *   <li>Ensures symmetry in the distance matrix and sets diagonal values to zero.</li>
     * </ul>
     *
     * @param initialProfiles a {@link LinkedList} of {@link Profile} objects representing the initial profiles
     *                        to be used in the algorithm. Each profile corresponds to a leaf node.
     */
    private void initialize(LinkedList<Profile> initialProfiles) {
        // put each sequence into a leaf node
        nodesOnMatrix = new Node[initialProfiles.size()];
        for (int i = 0; i < nodesOnMatrix.length; i++) {
            Node node = new Node(initialProfiles.get(i));
            nodesOnMatrix[i] = node;
        }

        // instantiate matrix
        distanceMatrix = new int[nodesOnMatrix.length][nodesOnMatrix.length];

        // for each Leaf, compute the distance to each other leaf
        for (int i = 0; i < nodesOnMatrix.length; i++) {
            for (int j = 0; j < nodesOnMatrix.length; j++) {
                if (j < i) {
                    // the d matrix is symmetric, so here we just place the already computed score into
                    // the correct index
                    distanceMatrix[i][j] = distanceMatrix[j][i];

                } else if (i == j) {
                    // diagonal is 0
                    distanceMatrix[i][j] = 0;

                }  else {
                    // for each combination of the initial sequences, compute their distance via their alignment score
                    Node node1 = nodesOnMatrix[i];
                    Node node2 = nodesOnMatrix[j];
                    String sequence1 = node1.getProfile().getInitialSequence();
                    String sequence2 = node2.getProfile().getInitialSequence();
                    int score = SequenceAlignment.computeAlignmentScore(sequence1, sequence2);
                    distanceMatrix[i][j] = score;
                }
            }
        }
    }

    /**
     * Computes the Neighbour Matrix for the current state of the distance matrix.
     *
     * <p>This method calculates the neighbour-joining values for each pair of nodes
     * in the distance matrix, which are used to determine the optimal nodes to merge
     * during the next step of the algorithm. The neighbour values are computed based
     * on the formula:
     * <pre>
     * N(i, j) = D(i, j) - (r(i) + r(j))
     * </pre>
     * where:
     * <ul>
     *   <li>D(i, j) is the distance between nodes i and j.</li>
     *   <li>r(i) is the average distance of node i to all other nodes.</li>
     *   <li>r(j) is the average distance of node j to all other nodes.</li>
     * </ul>
     *
     * <p>The neighbour matrix is symmetric and has a diagonal of zeros.
     * It is updated in-place during the computation.
     */
    private void computeNeighbourMatrix() {

        neighbourMatrix = new int[nodesOnMatrix.length][nodesOnMatrix.length];

        // for each Leaf, compute the distance to each other leaf
        for (int i = 0; i < nodesOnMatrix.length; i++) {
            for (int j = 0; j < nodesOnMatrix.length; j++) {

                if (j < i) {
                    // the d matrix is symmetric, so here we just place the already computed score into
                    // the correct index
                    neighbourMatrix[i][j] = neighbourMatrix[j][i];

                } else if (j == i) {
                    neighbourMatrix[i][j] = 0;

                }  else {
                    Node node1 = nodesOnMatrix[i];
                    Node node2 = nodesOnMatrix[j];
                    neighbourMatrix[i][j] = calculateNValue(node1, node2);
                }
            }
        }
    }

    /**
     * Calculates the neighbour-joining value (N-value) between two nodes.
     *
     * <p>The N-value is computed using the formula:
     * <pre>
     * N(i, j) = D(i, j) - (r(i) + r(j))
     * </pre>
     * where:
     * <ul>
     *   <li>N(i, j) is the neighbour-joining value between nodes i and j.</li>
     *   <li>D(i, j) is the distance between nodes i and j.</li>
     *   <li>r(i) is the average distance of node i to all other nodes.</li>
     *   <li>r(j) is the average distance of node j to all other nodes.</li>
     * </ul>
     *
     * <p>This value helps determine which pair of nodes should be merged next in the algorithm.
     *
     * @param node1 the first {@link Node}.
     * @param node2 the second {@link Node}.
     * @return the neighbour-joining value between the two nodes.
     */
    private int calculateNValue(Node node1, Node node2) {
        int rNode1 = computeDistanceMeanToEveryOtherNode(node1);
        int rNode2 = computeDistanceMeanToEveryOtherNode(node2);
        int distance = getDistanceBetween(node1, node2);

        return distance - (rNode1 + rNode2);
    }

    /**
     * Computes the average distance of the given node to all other nodes in the distance matrix.
     *
     * <p>This method calculates the mean of all distances between the specified node and
     * every other node, excluding the node itself. The formula used is:
     * <pre>
     * r(i) = sum(D(i, j)) / (N - 2)
     * </pre>
     * where:
     * <ul>
     *   <li>r(i) is the mean distance of node i to all other nodes.</li>
     *   <li>D(i, j) is the distance between nodes i and j.</li>
     *   <li>N is the total number of nodes.</li>
     * </ul>
     *
     * @param node the {@link Node} for which the mean distance is to be calculated.
     * @return the mean distance of the specified node to all other nodes.
     */
    private int computeDistanceMeanToEveryOtherNode(Node node) {
        int[] distances = getDistancesOf(node);
        int sum = 0;
        for (Integer distance : distances) {
            sum += distance;
        }
        return sum / (distanceMatrix.length - 2);
    }

    /**
     * Finds the pair of nodes with the smallest neighbour-joining value and combines them
     * into a new node.
     *
     * <p>This method iterates through the {@code neighbourMatrix} to identify the pair of nodes
     * with the minimum neighbour value, which indicates they are the closest in the context
     * of the neighbour-joining algorithm. The selected nodes are then merged into a new
     * {@link Node} with the selected nodes as its children.
     *
     * <p>The resulting node represents the combined structure of the two input nodes and will
     * be used in the next iteration of the algorithm.
     *
     * @return a new {@link Node} that has the two nearest nodes as its children.
     */
    private Node findAndCombineNearestNodes() {
        Node node1OfSmallestPair = null;
        Node node2OfSmallestPair = null;
        int smallestDistance = Integer.MAX_VALUE;
        for (int i = 0; i < nodesOnMatrix.length; i++) {
            for (int j = i; j < nodesOnMatrix.length; j++) {
                if (neighbourMatrix[i][j] < smallestDistance) {
                    smallestDistance = neighbourMatrix[i][j];
                    node1OfSmallestPair = nodesOnMatrix[i];
                    node2OfSmallestPair = nodesOnMatrix[j];
                }
            }
        }
        return new Node(node1OfSmallestPair, node2OfSmallestPair);
    }

    /**
     * Updates the distance matrix to include the new merged node and recalculates distances.
     *
     * <p>This method performs the following:
     * <ul>
     *   <li>Creates a new node array that excludes the two nodes that were merged and includes the new merged node.</li>
     *   <li>Initializes a new distance matrix to store the updated pairwise distances.</li>
     *   <li>Copies the distances between old nodes directly from the existing matrix, ensuring symmetry.</li>
     *   <li>Calculates the distances between the new merged node and the remaining nodes using the formula:
     *       <pre>
     *       D(new, k) = (D(i, k) + D(j, k) - D(i, j)) / 2
     *       </pre>
     *       where i and j are the merged nodes, and k is any other node.</li>
     * </ul>
     *
     * @param newNode the new {@link Node} that represents the merged structure of two existing nodes.
     */
    private void updateDistanceMatrix(Node newNode) {
        Node[] newNodes = new Node[nodesOnMatrix.length - 1];
        int[][] newDistanceMatrix = new int[newNodes.length][newNodes.length];

        Node mergedNode1 = newNode.getChildNode1();
        Node mergedNode2 = newNode.getChildNode2();

        int index = 0;
        for (Node oldNode : nodesOnMatrix) {
            if (!oldNode.equals(mergedNode1) && !oldNode.equals(mergedNode2)) {
                newNodes[index++] = oldNode;
            }
        }

        newNodes[newNodes.length - 1] = newNode;

        for (int i = 0; i < newNodes.length; i++) {
            for (int j = 0; j < newNodes.length; j++) {
                if (j < i) {
                    // symmetric filling
                    newDistanceMatrix[i][j] = newDistanceMatrix[j][i];

                } else if (i == j) {
                    newDistanceMatrix[i][j] = 0;

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
        nodesOnMatrix = newNodes;
        distanceMatrix = newDistanceMatrix;
    }

    /**
     * Computes the distance between an existing node and the newly merged node.
     *
     * <p>This method calculates the distance between an "old" node and a newly created
     * merged node based on the distances to the two child nodes of the merged node.
     * The formula used is:
     * <pre>
     * D(new, k) = (D(i, k) + D(j, k) - D(i, j)) / 2
     * </pre>
     * where:
     * <ul>
     *   <li>D(new, k) is the distance between the merged node and the existing node.</li>
     *   <li>D(i, k) and D(j, k) are the distances between the existing node and the two child nodes of the merged node.</li>
     *   <li>D(i, j) is the distance between the two child nodes of the merged node.</li>
     * </ul>
     *
     * @param node1 an existing {@link Node} for which the distance to the new node is to be calculated.
     * @param newNode the newly merged {@link Node} whose distance to the existing node is to be calculated.
     * @return the computed distance between the existing node and the new node.
     */
    private int computeNewDistanceBetween(Node node1, Node newNode) {
        Node node2 = newNode.getChildNode1();
        Node node3 = newNode.getChildNode2();

        return (this.getDistanceBetween(node1, node2) + this.getDistanceBetween(node1, node3) - this.getDistanceBetween(node2, node3)) / 2;
    }

    /**
     * @return distance between given two nodes in the distanceMatrix
     */
    private int getDistanceBetween(Node node1, Node node2) {
        int indexNode1 = 0;
        int indexNode2 = 0;

        for (int i = 0; i < nodesOnMatrix.length; i++) {
            if (nodesOnMatrix[i].equals(node1)) indexNode1 = i;
            if (nodesOnMatrix[i].equals(node2)) indexNode2 = i;
        }

        return distanceMatrix[indexNode1][indexNode2];
    }

    /**
     * @return the distances to every other node from the given node in the distanceMatrix
     */
    private int[] getDistancesOf(Node node) {
        int index = getIndexOfNode(node);
        return distanceMatrix[index];
    }

    /**
     * @return index of the given node in the nodeOnMatrix array
     */
    private int getIndexOfNode(Node node) {
        int index = 0;
        for (int i = 0; i < nodesOnMatrix.length; i++) {
            if (nodesOnMatrix[i].equals(node)) index = i;
        }
        return index;
    }

    /**
     * Represents a node in the phylogenetic tree constructed by the Neighbour Joining algorithm.
     *
     * <p>A node can either be:
     * <ul>
     *   <li>A leaf node, which corresponds to an initial profile with a specific sequence.</li>
     *   <li>An internal node, which represents the merging of two child nodes during the algorithm.</li>
     * </ul>
     *
     * <p>Each node maintains references to its parent node and its two child nodes (if applicable),
     * as well as associated profile information. Leaf nodes are directly associated with a profile,
     * while internal nodes represent merged profiles.
     *
     * <p>The key attributes of a node are:
     * <ul>
     *   <li>Parent node: The node that this node is connected to in the tree.</li>
     *   <li>Child nodes: The two nodes that were merged to form this node (if not a leaf).</li>
     *   <li>Profile: The profile or sequence data associated with this node.</li>
     *   <li>Name: A unique identifier or label for the node.</li>
     *   <li>Leaf flag: A boolean indicating whether the node is a leaf.</li>
     * </ul>
     *
     * <p>This class supports methods to retrieve and update its relationships and associated data.
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

        public boolean hasProfile() {
            return profile != null;
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
