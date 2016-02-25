package varun;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import com.sun.tools.javac.api.Formattable.LocalizedString;

import ij.ImageJ;
import net.imglib2.PointSampleList;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTree.KDTreeCursor;
import net.imglib2.KDTreeNode;
import net.imglib2.Localizable;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPointSampleList;
import net.imglib2.Sampler;
import net.imglib2.algorithm.kdtree.KDTreeNodeIterable;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.Cursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.MyKDtree.EucledianDistance;
import varun.MyKDtreeint.Distance;
import varun.MyKDtreeint.Node;

public class MyKDtree {

	/********* For an input image returns a PointSampleList ********/
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

	public static class Node<T> {

		private final int n;

		private final double medianValue;

		private final double[] nodePoint;

		private final int direction;

		private final PointSampleList<T> LeftTree;

		private final PointSampleList<T> RightTree;

		private final ArrayList<Point> newleftlist;
		
		private final ArrayList<Point> newrightlist;
		

		public Node(final double medianValue, double[] nodePoint, final int direction,
				final PointSampleList<T> LeftTree, final PointSampleList<T> RightTree,
				final ArrayList<Point> newleftlist,  ArrayList<Point> newrightlist
				) {
			assert LeftTree.numDimensions() == RightTree.numDimensions();
			this.n = LeftTree.numDimensions();

			this.nodePoint = nodePoint;

			this.medianValue = medianValue;
			this.direction = direction;
			this.LeftTree = LeftTree;
			this.RightTree = RightTree;

			this.newleftlist = newleftlist;
			
			this.newrightlist = newrightlist;
		

		}

		public int getnumDimensions() {
			return n;
		}

		public double getMedianValue() {
			return medianValue;

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

		public ArrayList<Point> getNewleftXlist() {
			return newleftlist;
		}


		public ArrayList<Point> getNewrightXlist() {
			return newrightlist;
		}

		

	}

	/********
	 * Constructor for the object Node that contains the Value at which a list
	 * is split up, the two split lists and the direction of the split
	 *********/

	public static class searchNode<T> {

		private final int n;

		public final double[] Position;

		protected PointSampleList<T> list;

		protected PointSampleList<T> finalsearchbranch;

		protected PointSampleList<T> finalnonsearchbranch;

		protected ArrayList<Point> Xlist;

		protected ArrayList<Point> Ylist;

		protected Node<T> finalnode;

		protected double Bestdistsquared;

		protected Node<T> farfinalnode;

		protected double Bestfardistsquared;

		public searchNode(final PointSampleList<T> list, final ArrayList<Point> Xlist, final ArrayList<Point> Ylist) {

			n = list.numDimensions();
			Position = new double[n];
			this.Xlist = Xlist;
			this.Ylist = Ylist;
			this.list = list;

		}

		public int getNumdimensions() {
			return n;
		}

		public void search(final RealLocalizable cursor, int direction) {
			cursor.localize(Position);
			Bestdistsquared = Double.MAX_VALUE;
			closestNode(list, Xlist, Ylist, direction);
		}

		public double getBestdist() {
			return Bestdistsquared;
		}

		public Node<T> getfinalnode() {
			return finalnode;
		}

		private void closestNode(final PointSampleList<T> list, ArrayList<Point> Xlist, ArrayList<Point> Ylist,
				int direction) {

			Node<T> currentBest = makeNode(list, Xlist, Ylist, direction);

			// System.out.println(" Closest Node: " + currentBest.medianValue);

			int nodedirection = currentBest.direction;

			double dist = 0;
			for (int d = 0; d < n; ++d) {

				dist += (Position[d] - currentBest.nodePoint[d]) * (Position[d] - currentBest.nodePoint[d]);
			}

			final double locationdiff = Position[currentBest.direction] - currentBest.nodePoint[currentBest.direction];
			final double axisdiff = locationdiff * locationdiff;
			final boolean leftbranchsearch = locationdiff < 0;

			final PointSampleList<T> searchBranch = leftbranchsearch ? currentBest.LeftTree : currentBest.RightTree;
			
			final ArrayList<Point> newXlist, newYlist, newnonsXlist, newnonsYlist;
			
			
			if (nodedirection == 0){
			 newXlist = leftbranchsearch ? currentBest.newleftlist : currentBest.newrightlist;
			newYlist = Ylist;
			
			newnonsXlist = leftbranchsearch ? currentBest.newleftlist : currentBest.newrightlist;
			 newnonsYlist = Ylist;
			
			}
			else{
				
			 newYlist = leftbranchsearch ? currentBest.newleftlist : currentBest.newrightlist;
			newXlist = Xlist;
			
			 newnonsYlist = leftbranchsearch ? currentBest.newleftlist : currentBest.newrightlist;
			 newnonsXlist = Xlist;
			
			}

			final PointSampleList<T> nonsearchBranch = leftbranchsearch ? currentBest.RightTree : currentBest.LeftTree;


			if (dist < Bestdistsquared) {
				Bestdistsquared = dist;
				finalnode = currentBest;

			}
			
			final boolean directionchoice = direction == n-1;
			final int otherdirection = directionchoice ? 0:direction +1;

			if ((searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1) > 2) {
				
				
				closestNode(searchBranch, newXlist, newYlist, otherdirection);
			}

			if (axisdiff <= Bestdistsquared
					&& (nonsearchBranch.realMax(otherdirection) - nonsearchBranch.realMin(otherdirection) + 1) > 2) {

			//	closestNode(nonsearchBranch, newnonsXlist, newnonsYlist, otherdirection);
			}
			
		}

		public double[] getMedian(ArrayList<Point> Xlist, ArrayList<Point> Ylist, int direction) {

			double[] medianPoint = new double[n];

			final boolean directionchoice = direction == n-1;
			
			
			
		final	ArrayList<Point> cordsort = directionchoice ?  Ylist:Xlist;

		final	ArrayList<Point> anticordsort = directionchoice ?  Ylist:Xlist;
	
		
		
		final int otherdirection = directionchoice ? 0:direction +1;
			
		int startindex = 0;
		
		int lastindex = cordsort.size()-1;
		
		int antilastindex = anticordsort.size()-1;
		
			

			
			
			

			medianPoint[direction] =cordsort.get(startindex).getDoublePosition(direction) + (cordsort.get(lastindex).getDoublePosition(direction) - cordsort.get(startindex).getDoublePosition(direction) )/2;

			
			medianPoint[otherdirection] = cordsort.get(startindex).getDoublePosition(otherdirection);
		

			
			System.out.println("Median Point x-cord "+medianPoint[0]);
			System.out.println("Median Point y-cord  "+medianPoint[1]);
			
			
			return medianPoint;

		}

		public Node<T> makeNode(PointSampleList<T> sortedlist, ArrayList<Point> Xlist, ArrayList<Point> Ylist,
				int direction) {

			

			/****
			 * To ward against running over the dimensionality, creating some
			 * local restrictions on the global variable direction
			 ****/
			if (direction == sortedlist.numDimensions())
				direction = 0;
			if ((sortedlist.realMax(direction) - sortedlist.realMin(direction) + 1) <= 2)
				return null;

			else {

				final boolean directionchoice = direction ==n-1;
				
				final int otherdirection = directionchoice ? 0:direction +1;
				
				

				double pivotElement;

				double[] point = new double[n];

				point = getMedian(Xlist, Ylist, direction);

				pivotElement = point[direction];

				final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
				final PointSampleList<T> RightTree = new PointSampleList<T>(n);

				final ArrayList<Point> newleftlist = new ArrayList<Point>(Xlist.size()/2);
				final ArrayList<Point> newrightlist = new ArrayList<Point>(Xlist.size()/2 +Xlist.size() % 2);
				
				final Cursor<T> listCursor = sortedlist.localizingCursor();
				
				
				
				
				
				final	ArrayList<Point> XorYlist = directionchoice ?  Ylist:Xlist;
				
				
				
				int index = 0;
				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point cord = new Point(listCursor);

					Point newpoint = XorYlist.get(index);

					

					if (listCursor.getDoublePosition(direction) < pivotElement) {

						LeftTree.add(cord, listCursor.get());
						newleftlist.add(newpoint);
						
					}

					else {

						RightTree.add(cord, listCursor.get());
						newrightlist.add(newpoint);
						

					}

					index++;
				}
				
				
				
				Node<T> node = new Node<T>(pivotElement, point, direction, LeftTree, RightTree, newleftlist,
						 newrightlist);

				return node;

			}

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

	/******
	 * Returns an arrayList of Points sorted in the given direction, useful for
	 * computing median
	 ******/

	public static <T extends RealType<T>> void sortpointList(ArrayList<Point> pointlist, int direction) {
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

	public static <T extends RealType<T>> PointSampleList<T> combineTrees(Node<T> Tree) {

		assert Tree.LeftTree.numDimensions() == Tree.RightTree.numDimensions();
		int n = Tree.LeftTree.numDimensions();

		PointSampleList<T> singleTree = new PointSampleList<T>(n);
		Cursor<T> treecursorA = Tree.LeftTree.cursor();
		Cursor<T> treecursorB = Tree.RightTree.cursor();

		while (treecursorA.hasNext()) {

			treecursorA.fwd();
			Point treepointA = new Point(n);

			treepointA.setPosition(treecursorA);

			singleTree.add(treepointA, treecursorA.get());

		}
		while (treecursorB.hasNext()) {

			treecursorB.fwd();
			Point treepointB = new Point(n);

			treepointB.setPosition(treecursorB);

			singleTree.add(treepointB, treecursorB.get());

		}

		return singleTree;

	}

	/****
	 * This is more of a "Java" way of getting NN by using an NN object and
	 * using getters and setters to get the Best Node point
	 ****/
	public static <T extends RealType<T>> void ConcisedistanceTransform(PointSampleList<BitType> list,
			ArrayList<Point> Xlist, ArrayList<Point> Ylist, PointSampleList<BitType> listzerosorones,
			RandomAccessibleInterval<T> imgout, final Distance dist) throws FileNotFoundException {

		// PrintStream out = new PrintStream(new
		// FileOutputStream("conKDtreemindist.txt"));
		// System.setOut(out);

		final Cursor<BitType> zerooronelistcursor = listzerosorones.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		final searchNode<BitType> Bestnode = new searchNode<BitType>(list, Xlist, Ylist);

		double distance = 0;

	//	while (zerooronelistcursor.hasNext()) {

			zerooronelistcursor.fwd();
			double mindistance = Double.MAX_VALUE;

			outbound.setPosition(zerooronelistcursor);

			Bestnode.search(zerooronelistcursor, 0);

			PointSampleList<BitType> singletree = combineTrees(Bestnode.getfinalnode());

			Cursor<BitType> singlecursor = singletree.cursor();

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				distance = dist.getDistance(zerooronelistcursor, singlecursor);
				mindistance = Math.min(distance, mindistance);

			}

			System.out.println(" Searchin points number: "+singletree.size());

			outbound.get().setReal(mindistance);

	//	}

		final Cursor<BitType> listcursor = list.localizingCursor();
		while (listcursor.hasNext()) {
			listcursor.fwd();
			outbound.setPosition(listcursor);
			outbound.get().setReal(0);

		}

	}
	
	
	

	public interface Distance {

		double getDistance(RealLocalizable cursor1, RealLocalizable cursor2);

		<T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint);

	}

	public static class EucledianDistance implements Distance {
		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {

			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);

			}

			return Math.sqrt(distance);

		}

		public <T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint) {
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

		public <T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += Math.abs(testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	/************
	 * Creating a bitType image from an image of type T by doing thresholding
	 * (not used currently)
	 ************/

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

	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());
		final Img<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());

		int n = bitimg.numDimensions();

		FloatType val = new FloatType(200);

		// ImageJFunctions.show(img).setTitle("KD-Tree input");

		createBitimage(img, bitimg, val);

		PointSampleList<FloatType> testlist = new PointSampleList<FloatType>(img.numDimensions());

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		IterableInterval<BitType> view = Views.interval(bitimg, new long[] { 0, 0 }, new long[] { 10, 10 });

		list = getList(bitimg);

		PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(n);

		PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(n);

		listonlyones = getvalueList(list, 1);
		listonlyzeros = getvalueList(list, 0);

		ArrayList<Point> Xpointsort = new ArrayList<Point>((int) listonlyones.size());
		ArrayList<Point> Ypointsort = new ArrayList<Point>((int) listonlyones.size());

		Xpointsort = getpointList(listonlyones);

		Ypointsort = getpointList(listonlyones);

		sortpointList(Xpointsort, 0); // Type points, sorted by X-coordinate
		sortpointList(Ypointsort, 1); // Type points, sorted by Y-coordinate

		// System.out.println(Xpointsort);
		// System.out.println(Ypointsort);

		long startTime = System.currentTimeMillis();
		ConcisedistanceTransform(listonlyones, Xpointsort, Ypointsort, listonlyzeros, imgout, new EucledianDistance());

		
		
		
		
		
		
		// new ImageJ();

//	ImageJFunctions.show(imgout).setTitle("KD-Tree output");
		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		//System.out.println(totalTime);

	}
}