package cmsc420_s23;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * This is a program that is a variation of the KDtree data structure. 
 * It has an internal polymorphic node class that has an internal and external node
 * Almost all of the methods that belong to the SMkdTree are implement using the 
 * internal and external classes. 
 * 
 * */
public class SMkdTree<LPoint extends LabeledPoint2D> {
	private Node root;
    private int rebuildOffset;
    private int deleteCount;
    private Rectangle2D rootCell;
    /**
     * This constructs a new (empty)
     * SMkdTree with the given rebuild offset and bounding box.
     * 
     * @param int rebuildOffset, Rectangle2D rootCell
     * */
	public SMkdTree(int rebuildOffset, Rectangle2D rootCell) { 
		this.rootCell = rootCell;
        this.rebuildOffset = rebuildOffset;
        root = new ExternalNode(0, null);
        root.cell = rootCell;
        deleteCount = 0;
	}
	
	/**
	 * This instantiates the abstract node class to be used for the InternalNode
	 * and ExternalNode classes 
	 * 
	 * */
	private abstract class Node {
		abstract LPoint find(Point2D q);
        abstract Node insert(LPoint pt, Rectangle2D cell) throws Exception;
        abstract Node delete(Point2D pt) throws Exception;
        abstract Node rebuild(Node node, Rectangle2D cell);
        abstract Node rebuildAfterInsert(Node node, Rectangle2D cell);
        abstract void buildPointList(List<LPoint> list);
        abstract LPoint nearestNeighbor(Point2D center, Rectangle2D rec, LPoint Best);
        abstract LPoint nearestNeighborVisit(Point2D center, Rectangle2D cell, LPoint best ,ArrayList<LPoint> visited);
        abstract ArrayList<String> list();
        abstract public String toString(); 
        
        Rectangle2D cell;
        int size;
        int count;
        int cutDim;
    }
	
	/**
	 * Defines the InternalNode class
	 * */
	public class InternalNode extends Node {
	    private double cutVal;
	    private Node left;//new ExternalNode(0, null);
	    private Node right;   //new ExternalNode(0, null);
	    
	    InternalNode(int cutDim, double cutVal, Node left, Node right) {
	        this.cutDim = cutDim; //Either zero which is a vertical cut or one if it is a horizontal cut
	        this.cutVal = cutVal; //The value at which the line is inserted for the cut. 
	        this.left = left;
	        this.right = right;
	        this.size = left.size + right.size; 
	    }
	    
	    /**
	     * This defines what to do for the find method when it is navigating internal nodes. Determines
	     * whether the left or right tree should be searched
	     * 
	     * @param Point2D q, q is the point being searched for
	     * @return LPoint of the q, returns null if q equals null
	     * */
	    LPoint find(Point2D q) {
	    	if (q == null) {
	    		return null;
	    	}
	        if (q.get(cutDim) < cutVal) {
	        	return left.find(q);
	        } 
	        else {
		        return right.find(q);
	        }
	    }
	    
	    /**
	     * This defines what the insert method should do when handling internal nodes. Determines if the
	     * LPoint pt should be inserted into the left or right tree. Once it goes down a side of the tree, 
	     * that trees size and insertion count are increased by one. 
	     * 
	     * @param LPoint pt, Rectangle2D cell 
	     * @return Node
	     * */
	    Node insert(LPoint pt, Rectangle2D cell) throws Exception {
	    	if (pt.get(cutDim) < cutVal) {
	    		left = left.insert(pt, cell.leftPart(cutDim, cutVal));
	        } 
	    	else {
	    		right = right.insert(pt, cell.rightPart(cutDim, cutVal));
	        }
	        size++;
	        count++;
	        return this;
	    }
	    
	    /**
	     * This function defines what the delete function should do when handling internal nodes. Determines
	     * if the left or right tree should be searched for the point to be deleted. At each subtree that is 
	     * searched their size is decremented
	     * 
	     * @param LPoint pt
	     * @return Node
	     * */
	    public Node delete(Point2D pt) throws Exception {
	    	if (pt.get(cutDim) < cutVal) {
	    		left = left.delete(pt);
	        } 
	    	else {
	    		right = right.delete(pt);
	        }
	    	size--;
	        return this;
	    } 
	    
	    /**
	     * When the insertion or deletion condition is met for a rebuild, this method takes care of gathering
	     * all of the points in that subtree and calls bulkCreate to rebuild the subtree. 
	     * 
	     * @param Node ndoe, Rectangle2D cell
	     * @return Reference to the new subtrees starting Node
	     * */
	    Node rebuild(Node node, Rectangle2D cell) {
	    	List<LPoint> points = new ArrayList<LPoint>();
	    	buildPointList(points);
	    	Node result = bulkCreate(points, cell);
	    	return result;
	    }
	    
	    /**
	     * After an insertion this method checks all of the subtrees to see if they should be rebuilt. 
	     * 
	     * @param Node node, Rectangle2D cell
	     * @return reference to the current node
	     * */
		Node rebuildAfterInsert(Node node, Rectangle2D cell) {
			if (node.count > ((node.size + rebuildOffset) / 2)) {
				return rebuild(node, cell);
			}
			else {
				left = left.rebuildAfterInsert(left, cell.leftPart(cutDim, cutVal));
				right = right.rebuildAfterInsert(right, cell.rightPart(cutDim, cutVal));
			}		
			return this;			
		}
		
		/**
		 * Recursively searches the left and right subtrees to find the external nodes in the subtree
		 * for the rebuild method. 
		 * */
		void buildPointList(List<LPoint> points) {
			left.buildPointList(points);
			right.buildPointList(points);
			
		} 
		
		/**
		 * This method gathers all of the points string outputs in right-to-left preorder traversal. 
		 * */
		ArrayList<String> list() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString()); // add this node
			list.addAll(right.list()); // add right
			list.addAll(left.list()); // add left
			return list;
		}
		/**
		 * Creates the toString for internal nodes
		 * */
		public String toString() {
			String prefix = (cutDim == 1) ? "y" : "x";
			return "(" + prefix + "=" + cutVal + ") " + size + ":" + count; 
		}

		/**
		 * The internal nodes nearest node determines which trees should and should not be checked
		 * by checking the distance from the center point to the cells of the other points. 
		 * 
		 * @param Point2D center, LPoint best
		 * @return LPoint closest point to center
		 * */
		LPoint nearestNeighbor(Point2D center, Rectangle2D cell, LPoint best) {
			Rectangle2D leftCell = cell.leftPart(cutDim, cutVal); // left childs cell
			Rectangle2D rightCell = cell.rightPart(cutDim, cutVal); // right childs cell
			
			if (center.get(cutDim) < cutVal) { // center is closer to left
				best = left.nearestNeighbor(center, leftCell, best); // search left subtree
				if (best != null) {
					if (rightCell.distanceSq(center) <= center.distanceSq(best.getPoint2D())) {// is right viable?
						best = right.nearestNeighbor(center, rightCell, best);
					}	
				}	
				else {
					if (rightCell.distanceSq(center) <= center.distanceSq(null)) {// is right viable?
						best = right.nearestNeighbor(center, rightCell, best);
					}	
				}		
			} 
			else { // center is closer to right
				best = right.nearestNeighbor(center, rightCell, best); // search right subtree
				if (best != null) {
					if (leftCell.distanceSq(center) <= center.distanceSq(best.getPoint2D())) {// is left viable?
						best = left.nearestNeighbor(center, leftCell, best);
					}			
				}
				else {
					if (leftCell.distanceSq(center) <= center.distanceSq(null)) {// is left viable?
						best = left.nearestNeighbor(center, leftCell, best);
					}
				}
			}
			return best;
		}

		/**
		 * This function determines which subtrees should be checked for the best cell. The visited parameter
		 * has added so this function in the external node can keep track of what points have been visited. 
		 * */
		LPoint nearestNeighborVisit(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited) {
			Rectangle2D leftCell = cell.leftPart(cutDim, cutVal); // left childs cell
			Rectangle2D rightCell = cell.rightPart(cutDim, cutVal); // right childs cell
			
			if (center.get(cutDim) < cutVal) { // center is closer to left
				best = left.nearestNeighborVisit(center, leftCell, best, visited); // search left subtree
				if (best != null) {
					if (rightCell.distanceSq(center) <= center.distanceSq(best.getPoint2D())) {// is right viable?
						best = right.nearestNeighborVisit(center, rightCell, best, visited);
					}	
				}	
				else {
					if (rightCell.distanceSq(center) <= center.distanceSq(null)) {// is right viable?
						best = right.nearestNeighborVisit(center, rightCell, best, visited);
					}	
				}
			} 
			else { // center is closer to right
				best = right.nearestNeighborVisit(center, rightCell, best, visited); // search right subtree
				if (best != null) {
					if (leftCell.distanceSq(center) <= center.distanceSq(best.getPoint2D())) {// is left viable?
						best = left.nearestNeighborVisit(center, leftCell, best, visited);
					}			
				}
				else {
					if (leftCell.distanceSq(center) <= center.distanceSq(null)) {// is left viable?
						best = left.nearestNeighborVisit(center, leftCell, best, visited);
					}
				}
			}
			return best;
		}
	}

	/**
	 * Defines the External Node class
	 * */
    public class ExternalNode extends Node {
        private LPoint point;
        
        /**
         * Constructor for the External Node class
         * */
        ExternalNode(int cutDim, LPoint point) {
            this.cutDim = cutDim;
            this.point = point;
            if (point != null) {
                size = 1;
            }
            else {
            	size = 0;
            }
        }
        
        /**
         * If the current node is an external node it is compared to see if it null or 
         * matches q. if it does not match q it returns null. 
         * */
        LPoint find(Point2D q) {
        	if (point == null) {
        		return null;
        	}
            if (point.getX() == q.getX() && point.getY() == q.getY()) {
                return point;
            } 
            else {
                return null;
            }
        }
        
        /**
         * The insert function for an external node checks to see if the current point is null, if 
         * it is then bulkCreate is called just for point pt. If the point is not null then both of 
         * them are put into bulkCreate so they be added in the correct order.	
         * 
         * @param LPoint pt, Rectangle2D cell
         * @exception if pt already exists in the tree an exception is thrown
         * @return reference to the starting node of the new node or subtree
         * */
        public Node insert(LPoint pt, Rectangle2D cell) throws Exception { 
            List<LPoint> points = new ArrayList<LPoint>();
            if (point == null) {
        		points.add(pt);
                return bulkCreate(points, cell);
            }
        	else {
            	if (point.getPoint2D().equals(pt.getPoint2D())) { // if the point already exists
                    throw new Exception("Insertion of duplicate point");
                }    
                // create a list of the two points and rebuild the tree
                points.add(pt);
                points.add(point);
                return bulkCreate(points, cell);
        	}
           
        }
        
        /**
         * If the pt matches point then the point is deleted by setting it to null and the deletion counter
         * is incremented by one. 
         * 
         * @param Point2D pt
         * @exception if the point is null or does not match pt then an error is thrown 
         * @return reference to the current node
         * */
        Node delete(Point2D pt) throws Exception {
            if (point == null || !pt.equals(point.getPoint2D())) {
				throw new Exception("Deletion of nonexistent point");
			}
            point = null;
	        deleteCount++;
            return this;
        }
        
        /**
         * This function does nothing for external nodes
         * */
        Node rebuild(Node node, Rectangle2D cell) {
        	return this;
        }
        
        /**
         * This function does nothing for external nodes
         * */
		Node rebuildAfterInsert(Node node, Rectangle2D cell) {
        	return this;
		}

		/**
		 * If point is non null then it is added to the list of points
		 * */
		void buildPointList(List<LPoint> points) {
			if (point != null) {
				points.add(point);
			}
		} 
		
		/**
		 * Adds the current point toString to the list 
		 * */
		ArrayList<String> list() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString()); // add this node
			return list;
		}
		
		/**
		 * Returns the string for the current point
		 * */
		public String toString() {
			if (point == null) {
				return ("[null]");
			}
			return ("[" + point.toString() +"]"); 
		}
		
		/**
		 * If the point is null then best is returned. If best is null and point is not then best is set equal 
		 * to the current point. If neither applies then the current point is compared to best point to see which
		 * point is closer and adjusts the variable best accordingly. 
		 * 
		 * @return best 
		 * */
		LPoint nearestNeighbor(Point2D center, Rectangle2D cell, LPoint best) {
			if (point == null) {
				return best;
			}
			if (best == null) {
				best = point;
			}
			else if (center.distance(point.getPoint2D()) <= center.distance(best.getPoint2D())) {
				if (center.distance(point.getPoint2D()) == center.distance(best.getPoint2D())) {
					if (point.getX() < best.getX()) {
						best = point;
					}
				}
				else {
					best = point;	
				}
			}
			return best;
		}
		
		/**
		 * Functions the same as nearest neighbor determining the current best. At the end of the function the 
		 * current point is added to the visited list regardless of what happens to best.   
		 * */
		LPoint nearestNeighborVisit(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited) {
			if (point == null) {
				return best;
			}
			if (best == null) {
				best = point;
			}
			else if (center.distance(point.getPoint2D()) <= center.distance(best.getPoint2D())) {
				if (center.distance(point.getPoint2D()) == center.distance(best.getPoint2D())) {
					if (point.getX() < best.getX()) {
						best = point;
					}
				}
				else {
					best = point;	
				}
			}
			visited.add(point);
			return best;
		}	
    }
    
    /** 
     * This resets the tree to its initial (empty) condition.
     * */
	public void clear() {
		root = new ExternalNode(0, null);
	    root.size = 0;
	    deleteCount = 0; 
	}
	
	/** 
	 * Returns the number of points in the tree. 
	 * */
	public int size() { 
		 return root.size;
	}
	
	/**
	 * Returns the current value of the deletion counter, described above. 
	 * */
	public int deleteCount() { 
		 return deleteCount;
	}
	
	/**
	 * Calls the find function in the internal and external node. 
	 * */
	public LPoint find(Point2D q) {
	    return root.find(q);
	}
	
	/**
	 * Checks to see if pt is null. If not it calls the delete function in the internal and
	 * external node classes. It then checks the to see if the rebuild condition has been met 
	 * and if so the rebuild function is called and deletion count is set to 0.  
	 * */
	public void delete(Point2D pt) throws Exception { 
		if (pt == null) {
			return;
		}
		root.delete(pt);
		if (deleteCount > root.size) {
			root = root.rebuild(root, rootCell);
			deleteCount = 0;
		}
	}
    /**
     * Throws an exception if the point pt lies outside of the bounds of RootCell. 
     * Calls insert method and rebuildAfterInsert method for internal and external node. 
     * */
    public void insert(LPoint pt) throws Exception {
    	if (!rootCell.contains(pt.getPoint2D())) {
            throw new Exception("Attempt to insert a point outside bounding box");
        }
        root = root.insert(pt, rootCell);
        root = root.rebuildAfterInsert(root, rootCell);
    }
    
    /**
     * Calls the list function in internal and external node class. 
     * */
    public ArrayList<String> list() {
	    return root.list();
	}
    
    /**
     * Calls the nearestNeighbor function in the internal and external node class. 
     * */
	public LPoint nearestNeighbor(Point2D center) { 
		return root.nearestNeighbor(center, rootCell, null);
	}
	
	/**
	 * Calls the nearestNeighborVisit function in the internal and external node class. 
	 * Then it uses the comparator to sort the list. 
	 * */
	public ArrayList<LPoint> nearestNeighborVisit(Point2D center) { 
		ArrayList<LPoint> visited = new ArrayList<LPoint>();
		root.nearestNeighborVisit(center, rootCell, null, visited); 
		visited.sort((Comparator<? super LPoint>) new ByXThenY());
		return visited;
	}
    
	/**
	 * If there are no points in list than an empty ExternalNode is returned. If there is one point
	 * then an ExternalNode with that point as a parameter is returned. If there are two or more 
	 * points in the list pts, then the helpter method cutDimension is called to determine the 
	 * cutDimension for the list of points. It then sorts them based off the cutDimesion. The cut 
	 * value is then determined by the helper method cutValue for the list pts. Then a for loop is 
	 * called to determine how many points should be sorted into the left and right subtree. Then 
	 * bulkCreate is recursively called on the left and right subtrees to be created. 
	 * 
	 * */
    public Node bulkCreate(List<LPoint> pts, Rectangle2D cell) {
        int n = pts.size();
        if (n == 0) {
            return new ExternalNode(0, null);
        } 
        else if (n == 1) {
            return new ExternalNode(0, pts.get(0));
        }
        else {
        	int cutDim = cutDimension(pts, cell);
        	double cutVal;
        	
        	if (cutDim == 0) {
        		pts.sort((Comparator<? super LPoint>) new ByXThenY());
        	}
        	else {
        		pts.sort((Comparator<? super LPoint>) new ByYThenX());
        	}
        	
        	if (cutDim == 0) {
        		if (pts.get(0).getX() == (pts.get(n - 1).getX())) {
            		cutDim = 1;
            	}
        	}
        	else {
        		if (pts.get(0).getY() == (pts.get(n - 1).getY())) {
            		cutDim = 0;
            	}
        	}

        	cutVal = cutValue(cutDim, pts, cell);
        	int j = 0;
        	for (int i = 0; i < pts.size(); i++) {
        		if (cutDim == 0) {
        			if (pts.get(i).getX() >= cutVal) {
        				break;
        			}
        		}
        		else {
        			if (pts.get(i).getY() >= cutVal) {
        				break;
        			}
        		}
        		j++;
        	}
        	
        	Node left = bulkCreate(pts.subList(0, j), cell.leftPart(cutDim, cutVal));
        	Node right = bulkCreate(pts.subList(j, pts.size()), cell.rightPart(cutDim, cutVal));
        	
        	return new InternalNode (cutDim, cutVal, left, right);
        }
        
    }
    
    /**
     * Determines the cutting dimension by comparing the x and y sides to see which is longer. If
     * there are equal then the cutting dimension is 0. 
     * */
    public int cutDimension(List<LPoint> pts, Rectangle2D cell) {
    	if (cell.high.getX() - cell.low.getX() >= cell.high.getY() - cell.low.getY()) {
    		return 0;
    	}
    	return 1;
    }
    
    /**
     * This function determines the cutValue for the list of points by the cell parameter. The cutVal
     * is determined by finding the midpoint of the longer axis. If the points all lay on one side of 
     * the cutValue then the cutValue slides to closest point. 
     * */
    public double cutValue(int cutDim, List<LPoint> pts, Rectangle2D cell) {
    	double cutVal = 0;
    	if (cutDim == 0) {
    		cutVal = ((cell.high.getX() + cell.low.getX()) / 2);
    		if (pts.get(0).getX() > cutVal) {
    			return cutVal = pts.get(0).getX();
    		}
    		if (pts.get(pts.size() - 1).getX() < cutVal) {
    			return cutVal = pts.get(pts.size() - 1).getX();
    		}
    	}
    	else {
    		cutVal = ((cell.high.getY() + cell.low.getY()) / 2);	
    		if (pts.get(0).getY() > cutVal) {
    			return cutVal = pts.get(0).getY();
    		}
    		if (pts.get(pts.size() - 1).getY() < cutVal) {
    			return cutVal = pts.get(pts.size() - 1).getY();
    		}
    	}
    	return cutVal;
    }
    
    /**
     * Sorts the list by the X value first and then the Y value. 
     * */
    private class ByXThenY implements Comparator<LPoint> {
        public int compare(LPoint pt1, LPoint pt2) {
            if (pt1.getX() < pt2.getX()) {
                return -1;
            } else if (pt1.getX() > pt2.getX()) {
                return 1;
            } else {
                if (pt1.getY() < pt2.getY()) {
                    return -1;
                } else if (pt1.getY() > pt2.getY()) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }
    }

    /**
     * Sorts the list by the Y value and then X value
     * */
    private class ByYThenX implements Comparator<LPoint> {
        public int compare(LPoint pt1, LPoint pt2) {
            if (pt1.getY() < pt2.getY()) {
                return -1;
            } else if (pt1.getY() > pt2.getY()) {
                return 1;
            } else {
                if (pt1.getX() < pt2.getX()) {
                    return -1;
                } else if (pt1.getX() > pt2.getX()) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }
    }

	// The following is needed only for the Challenge Problem

	private class LPointIterator implements Iterator<LPoint> {
		public LPointIterator() { /* ... */ }
		public boolean hasNext() { /* ... */ return false; }
		public LPoint next() throws NoSuchElementException { /* ... */ return null; }

	}
	
	public LPointIterator iterator() { /* ... */ return new LPointIterator(); }
}

