
/**
 * A node in the BnB-tree. Each node corresponds to a path in the problem-graph. 
 * This path can be reconstructed by looking at the edge-field and then recursively 
 * on the edge-fields of the nodes-parent. The lower bound field is set in the 
 * BnB-class.
 * 
 */
public class BnBNode {
	final BnBNode parent;
	final Edge edge;
	double lowerBound;

	/** Construct the first node */
	BnBNode(Edge edge){
		this(null,edge);
	}
	/** Construct any of the subsequent nodes */
	BnBNode(BnBNode parent, Edge edge){
		this.parent = parent;
		this.edge = edge;
	}

	public String toString(){
		return String.format("BnBNode[%s,lowerBound=%.2f,id=%d]",edge.toString(),lowerBound);
	}
}
