package varun;

import java.awt.List;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import javax.swing.tree.TreeNode;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.MyKDtree.Distance;
import varun.MyKDtree.Node;

public class Tree {

	public static class Node<T> {
		private final int n;

		private final double[] nodePoint;

		private final int direction;

		public final PointSampleList<T> LeftTree;

		public final PointSampleList<T> RightTree;

		public Node(double[] nodePoint, final int direction, final PointSampleList<T> LeftTree,
				final PointSampleList<T> RightTree) {

			assert LeftTree.numDimensions() == RightTree.numDimensions();
			this.n = LeftTree.numDimensions();

			this.nodePoint = nodePoint;

			this.direction = direction;
			this.LeftTree = LeftTree;
			this.RightTree = RightTree;

		}

		public Node() {
			this.n = 0;
			this.nodePoint = null;
			this.direction = 0;
			this.LeftTree = null;
			this.RightTree = null;
			
			
		}

		public int getnumDimensions() {
			return n;
		}

		public double getMedianValue() {
			return nodePoint[direction];

		}

		public int getDirection() {
			return direction;
		}

		public PointSampleList<T> getLeftTree() {
			return LeftTree;
		}

		public PointSampleList<T> getRightTree() {
			return RightTree;
		}

		public void localize(double[] point) {

			for (int d = 0; d < n; ++d)
				point[d] = (double) nodePoint[d];
		}

	}

	public static class fullTree<T> {

		private final Node<T> parentnode;

		private final Node<T> leftnode;

		private final Node<T> rightnode;

		private final int treeindex;

		public fullTree(Node<T> parentnode, Node<T> leftnode, Node<T> rightnode, int treeindex) {

			this.parentnode = parentnode;
			this.leftnode = leftnode;
			this.rightnode = rightnode;
			this.treeindex = treeindex;

		}

		public int getTreeindex() {
			return treeindex;
		}

		public Node<T> getParentnode() {
			return parentnode;
		}

		public Node<T> getLeftnode() {
			return leftnode;
		}

		public Node<T> getRightnode() {
			return rightnode;
		}

	}

	/********
	 * Constructor for the object Node that contains the Value at which a list
	 * is split up, the two split lists and the direction of the split
	 *********/

	public static class searchNode<T> {

		private final int n;

		public final double[] Position;

		protected final ArrayList<Node<T>> allnodes;

		protected Node<T> finalnode;

		protected double Bestdistsquared;

		protected double finallocationdiff;

		protected double Bestaxisdiffsquared;

		public searchNode(final ArrayList<Node<T>> allnodes) {

			n = allnodes.get(0).getnumDimensions();
			Position = new double[n];

			this.allnodes = allnodes;

		}

		public int getNumdimensions() {
			return n;
		}

		public void search(final RealLocalizable cursor, final int direction) throws FileNotFoundException {
			cursor.localize(Position);
			Bestdistsquared = Double.MAX_VALUE;
			finallocationdiff = Double.MAX_VALUE;
			closestNode(allnodes.get(allnodes.size() - 1), direction);
		}

		public double getBestdist() {
			return Bestdistsquared;
		}

		public double getBestaxisdiffdist() {
			return Bestaxisdiffsquared;
		}

		public Node<T> getfinalnode() {
			return finalnode;
		}

		private void closestNode(final Node<T> tree, final int direction) throws FileNotFoundException {

			final boolean directionchoice = direction == n - 1;
			final int otherdirection = directionchoice ? 0 : direction + 1;

			final Node<T> currentBest = tree;

			double dist = 0;

			for (int d = 0; d < n; ++d) {

				dist += Math.pow((Position[d] - currentBest.nodePoint[d]), 2);
			}

			final double locationdiff = Position[currentBest.direction] - currentBest.nodePoint[currentBest.direction];

			double axisdiff = locationdiff * locationdiff;

			final boolean leftbranchsearch = locationdiff < 0;

			if (dist <= Bestdistsquared) {

				Bestdistsquared = dist;
				finalnode = currentBest;

			}

			final PointSampleList<T> searchBranch = leftbranchsearch ? currentBest.LeftTree : currentBest.RightTree;

			final PointSampleList<T> nonsearchBranch = leftbranchsearch ? currentBest.RightTree : currentBest.LeftTree;

			final ArrayList<Point> newXlist, newYlist, newnonsXlist, newnonsYlist;

			final boolean nodedirectionchoice = currentBest.direction == n - 1;

			// closestNode(searchBranch, newXlist, newYlist, otherdirection,
			// dist);

			// if ( axisdiff <= Bestdistsquared)

			// closestNode(nonsearchBranch, newnonsXlist, newnonsYlist,
			// otherdirection, dist);
		}

	}

	public static <T extends RealType<T>> void makeTree(final Node<T> mainnode,ArrayList<fullTree<T>> treelist,
			final int direction,final int treeindex) {

		

			int n = mainnode.getLeftTree().numDimensions();

			final boolean directionchoice = direction == n - 1;
			final int otherdirection = directionchoice ? 0 : direction + 1;

			

			 Node<T> childnodeA = new Node<T>();
			 Node<T> childnodeB = new Node<T>();	
					
					
					childnodeA = makeNode(mainnode.getLeftTree(), otherdirection);
			        childnodeB = makeNode(mainnode.getRightTree(), otherdirection);
			      
		final	fullTree<T> tree = new fullTree<T>(mainnode, childnodeA, childnodeB, treeindex);
		
		

		if(childnodeA!=null)
			makeTree(childnodeA,treelist,direction,treeindex+1);
		if(childnodeB!=null)
			makeTree(childnodeB, treelist, direction, treeindex+2);
		
		treelist.add(tree);
			
System.out.println("treeindex :"+tree.treeindex);


		
	}

	/*********
	 * Here I return an Arraylist of Node<T> type by a clever trick, I give in
	 * an ArrayList of Node<T> containing only the rootNode and in the course of
	 * the routine below overwrite that list with the nodes of all the subtrees.
	 * The list then contains the object Node<T> for all the subtrees including
	 * the rootTree. Is this clever or dangerous way to do it?
	 ********/
	public static <T extends RealType<T>> void getTree(PointSampleList<T> list, ArrayList<Node<T>> allnodes,
			int direction) {

		if (list.realMax(direction) - list.realMin(direction) + 1 <= 2)
			return;

		else {
			int n = list.numDimensions();

			final boolean directionchoice = direction == n - 1;
			final int otherdirection = directionchoice ? 0 : direction + 1;

			ArrayList<Point> Xpointsort = new ArrayList<Point>();
			ArrayList<Point> Ypointsort = new ArrayList<Point>();

			Xpointsort = getpointList(list);

			Ypointsort = getpointList(list);

			sortpointList(Xpointsort, 0); // Type points, sorted by X-coordinate
			sortpointList(Ypointsort, 1); // Type points, sorted by Y-coordinate

			final ArrayList<Point> cordsort = directionchoice ? Ypointsort : Xpointsort;

			final ArrayList<Point> anticordsort = directionchoice ? Xpointsort : Ypointsort;

			if (cordsort.size() > 2) {

				final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
				final PointSampleList<T> RightTree = new PointSampleList<T>(n);

				double[] point = new double[n];
				point = getMedian(cordsort, anticordsort, direction, n);

				// point = getMean(list, direction);

				final Cursor<T> listCursor = list.localizingCursor();

				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point cord = new Point(listCursor);

					if (listCursor.getDoublePosition(direction) < point[direction])

						LeftTree.add(cord, listCursor.get());

					else

						RightTree.add(cord, listCursor.get());

				}

				Node<T> node = new Node<T>(point, direction, LeftTree, RightTree);

				getTree(LeftTree, allnodes, otherdirection);

				getTree(RightTree, allnodes, otherdirection);

				nodetoList(node, allnodes);

			}

			else
				getTree(list, allnodes, otherdirection);

		}

	}

	public static <T extends RealType<T>> void nodetoList(final Node<T> node, final ArrayList<Node<T>> allnodes) {
		allnodes.add(node);
	}

	public static <T extends RealType<T>> double[] getMedian(ArrayList<Point> cordsort, ArrayList<Point> anticordsort,
			int direction, int n) {

		final boolean directionchoice = direction == n - 1;
		final int otherdirection = directionchoice ? 0 : direction + 1;

		final double[] medianPoint = new double[n];

		int medianindexA = (cordsort.size() - 1) / 2;

		medianPoint[direction] = (cordsort.get(medianindexA).getDoublePosition(direction));

		int medianindexB = (anticordsort.size() - 1) / 2;

		medianPoint[otherdirection] = (anticordsort.get(medianindexB).getDoublePosition(otherdirection));

		return medianPoint;

	}

	public static <T extends RealType<T>> double[] getMean(PointSampleList<T> list, int direction) {
		int n = list.numDimensions();

		final boolean directionchoice = direction == n - 1;
		final int otherdirection = directionchoice ? 0 : direction + 1;

		final double[] meanPoint = new double[n];

		int medianindexA = (int) (list.size() - 1) / 2;
		Cursor<T> listCursor = list.localizingCursor();
		listCursor.jumpFwd(medianindexA);

		meanPoint[direction] = (listCursor.getDoublePosition(direction));

		meanPoint[otherdirection] = (listCursor.getDoublePosition(otherdirection));

		return meanPoint;

	}

	public static <T extends RealType<T>> Node<T> makeNode(PointSampleList<T> list, int direction) {
		int n = list.numDimensions();

		final boolean directionchoice = direction == n - 1;
		final int otherdirection = directionchoice ? 0 : direction + 1;

		ArrayList<Point> Xpointsort = new ArrayList<Point>();
		ArrayList<Point> Ypointsort = new ArrayList<Point>();

		Xpointsort = getpointList(list);

		Ypointsort = getpointList(list);

		sortpointList(Xpointsort, 0); // Type points, sorted by X-coordinate
		sortpointList(Ypointsort, 1); // Type points, sorted by Y-coordinate

		final ArrayList<Point> cordsort = directionchoice ? Ypointsort : Xpointsort;

		final ArrayList<Point> anticordsort = directionchoice ? Xpointsort : Ypointsort;
		if (cordsort.size() > 2) {
		double[] point = new double[n];

		point = getMedian(cordsort, anticordsort, direction, n);

		final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
		final PointSampleList<T> RightTree = new PointSampleList<T>(n);

		final Cursor<T> listCursor = list.localizingCursor();

		while (listCursor.hasNext()) {

			listCursor.fwd();

			Point cord = new Point(listCursor);

			if (listCursor.getDoublePosition(direction) < point[direction])

				LeftTree.add(cord, listCursor.get());

			else

				RightTree.add(cord, listCursor.get());

		}

		Node<T> node = new Node<T>(point, direction, LeftTree, RightTree);

		return node;
		}
		else
			return null;
	}

	public static void sortpointList(ArrayList<Point> pointlist, int direction) {
		if (pointlist.size() <= 1)
			return;

		else {

			// the first element belonging to the right list childB
			final int splitIndex = (int) pointlist.size() / 2;

			Iterator<Point> iterator = pointlist.iterator();

			final ArrayList<Point> childA = new ArrayList<Point>((int) pointlist.size() / 2);

			final ArrayList<Point> childB = new ArrayList<Point>((int) (pointlist.size() / 2 + pointlist.size() % 2));

			int index = 0;

			while (iterator.hasNext()) {
				iterator.next();

				if (index < splitIndex)
					childA.add(pointlist.get(index));

				else

					childB.add(pointlist.get(index));

				index++;

			}

			sortpointList(childA, direction);

			sortpointList(childB, direction);

			mergepointListValue(pointlist, childA, childB, direction);

			/********
			 * The part below removes the duplicate entries in the sorted array
			 ********/
			int j = 0;

			for (int i = 0; i < pointlist.size(); ++i) {

				j = i + 1;
				while (j < pointlist.size()) {

					if (pointlist.get(i).getDoublePosition(direction) == pointlist.get(j)
							.getDoublePosition(direction)) {

						pointlist.remove(j);

					}

					else {
						++j;
					}

				}

			}

		}

	}

	/// ***** Returns a sorted list *********////
	public static void mergepointListValue(ArrayList<Point> sortedlist, ArrayList<Point> listA, ArrayList<Point> listB,
			int direction) {

		int i = 0, j = 0, k = 0;

		while (i < listA.size() && j < listB.size()) {

			if (listA.get(i).getDoublePosition(direction) < (listB.get(j).getDoublePosition(direction))) {

				sortedlist.set(k, listA.get(i));

				++i;
				++k;
			}

			else {

				sortedlist.set(k, listB.get(j));

				++j;
				++k;

			}

		}

		while (i < listA.size()) {
			sortedlist.set(k, listA.get(i));
			++i;
			++k;

		}

		while (j < listB.size()) {
			sortedlist.set(k, listB.get(j));
			++j;
			++k;

		}

	}

	public static <T extends RealType<T>> ArrayList<Point> getpointList(PointSampleList<T> sortedlist) {
		ArrayList<Point> pointlist = new ArrayList<Point>();

		Cursor<T> newlistCursor = sortedlist.localizingCursor();

		while (newlistCursor.hasNext()) {

			newlistCursor.fwd();

			Point newcord = new Point(newlistCursor);

			pointlist.add(newcord);

		}

		return pointlist;
	}

	public static <T extends RealType<T>> void createBitimage(RandomAccessibleInterval<T> img, Img<BitType> imgout,
			T ThresholdValue) {

		final Cursor<T> bound = Views.iterable(img).localizingCursor();

		final RandomAccess<BitType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().compareTo(ThresholdValue) > 0) {

				outbound.get().setOne();

			}

			else {

				outbound.get().setZero();

			}

		}
	}

	public static <T extends RealType<T>> PointSampleList<T> getList(IterableInterval<T> img) {

		int n = img.numDimensions();

		final Cursor<T> first = img.cursor();

		// A point sample list with coordinates declared and initialized.
		PointSampleList<T> parent = new PointSampleList<T>(n);

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(n);

			cord.setPosition(first);

			parent.add(cord, first.get().copy());

		}

		return parent;

	}

	public static PointSampleList<BitType> getvalueList(PointSampleList<BitType> list, int val) {

		int n = list.numDimensions();

		final Cursor<BitType> first = list.cursor();

		// A point sample list with coordinates declared and initialized.
		PointSampleList<BitType> parent = new PointSampleList<BitType>(n);

		while (first.hasNext()) {
			first.fwd();
			if (first.get().getInteger() == val) {
				Point cord = new Point(n);

				cord.setPosition(first);

				parent.add(cord, first.get().copy());

			}
		}
		return parent;

	}

	public interface Distance {

		double getDistance(RealLocalizable cursor1, RealLocalizable cursor2);

		double getDistance(RealLocalizable listcursor, double[] testpoint);

	}

	public static class EucledianDistance implements Distance {
		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {

			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);

			}

			return Math.sqrt(distance);

		}

		public double getDistance(RealLocalizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += (testpoint[d] - listcursor.getDoublePosition(d))
						* (testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	public static class MannhattanDistance implements Distance {

		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {
			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.abs(cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d));

			}

			return distance;
		}

		public double getDistance(RealLocalizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += Math.abs(testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());

		int n = bitimg.numDimensions();

		FloatType val = new FloatType(200);

		// ImageJFunctions.show(img).setTitle("KD-Tree input");

		createBitimage(img, bitimg, val);

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		IterableInterval<BitType> view = Views.interval(bitimg, new long[] { 0, 0 }, new long[] { 20, 20 });

		list = getList(bitimg);

		PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(n);

		PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(n);

		listonlyones = getvalueList(list, 1);
		listonlyzeros = getvalueList(list, 0);

		ArrayList<Node<BitType>> allnodes = new ArrayList<Node<BitType>>();

		// getTree(listonlyones, allnodes, 0);

		final Node<BitType> rootnode = makeNode(list, 0);

		final ArrayList<fullTree<BitType>> testtree = new ArrayList<fullTree<BitType>>();

		makeTree(rootnode,testtree, 0,0);
		
		System.out.println(testtree.size());
		
	//	for (int index = 0; index < testtree.size(); ++index)
		//	System.out.println(testtree.get(index).treeindex);

		/*
		 * 
		 * 
		 * PrintStream out = new PrintStream(new
		 * FileOutputStream("output.txt")); System.setOut(out);
		 * 
		 * Cursor<BitType> listcursor = listonlyones.localizingCursor();
		 * while(listcursor.hasNext()){ listcursor.fwd(); System.out.println(
		 * " List points X cordinates:"+listcursor.getDoublePosition(0));
		 * System.out.println(" List points Y coordinates:"
		 * +listcursor.getDoublePosition(1)); }
		 * 
		 * for (int index = allnodes.size()-1; index>=0; --index){
		 * System.out.println("Median root node: "
		 * +allnodes.get(index).getMedianValue()); System.out.println("Index: "
		 * +index); System.out.println("Direction: "
		 * +allnodes.get(index).getDirection());
		 * 
		 * Cursor<BitType> treecursorleft =
		 * allnodes.get(index).LeftTree.localizingCursor();
		 * 
		 * while(treecursorleft.hasNext()){ treecursorleft.fwd();
		 * 
		 * 
		 * System.out.println("Left Tree points X: "
		 * +treecursorleft.getDoublePosition(0)); System.out.println(
		 * "Left Tree points Y: "+treecursorleft.getDoublePosition(1));
		 * 
		 * 
		 * } }
		 */

	}

}
