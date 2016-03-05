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
import net.imglib2.KDTreeNode;
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
		
		private final int treeindex;

		public final PointSampleList<T> Leftsublist;

		public final PointSampleList<T> Rightsublist;
		
		public Node<T> leftchild;
		
		public  Node<T> rightchild;

		public Node(double[] nodePoint, final int direction, final int treeindex, final PointSampleList<T> Leftsublist,
				final PointSampleList<T> Rightsublist, final Node<T> leftchild, final Node<T> rightchild) {

			assert Leftsublist.numDimensions() == Rightsublist.numDimensions();
			this.n = Leftsublist.numDimensions();

			this.nodePoint = nodePoint;

			this.treeindex = treeindex;
			this.direction = direction;
			this.Leftsublist = Leftsublist;
			this.Rightsublist = Rightsublist;
			this.leftchild = leftchild;
			this.rightchild = rightchild;

		}

		public Node(final Node< T > node )
		{
			this.n = node.n;			
			this.nodePoint = node.nodePoint;
			this.direction = node.direction;
			this.Leftsublist = node.Leftsublist;
			this.Rightsublist = node.Rightsublist;
			this.leftchild = node.leftchild;
			this.rightchild = node.rightchild;
			this.treeindex = node.treeindex;
			
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

		public PointSampleList<T> getLeftsublist() {
			return Leftsublist;
		}

		public PointSampleList<T> getRightsublist() {
			return Rightsublist;
		}

		public void localize(double[] point) {

			for (int d = 0; d < n; ++d)
				point[d] = (double) nodePoint[d];
		}

		public int getTreeindex() {
			return treeindex;
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

	public static <T extends RealType<T>> Node<T> makeNode(PointSampleList<T> list, int direction, int treeindex) {
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
	
		double[] point = new double[n];

		point = getMedian(cordsort, anticordsort, direction, n);

		final PointSampleList<T> Leftsublist = new PointSampleList<T>(n);
		final PointSampleList<T> Rightsublist = new PointSampleList<T>(n);

		final Cursor<T> listCursor = list.localizingCursor();

		while (listCursor.hasNext()) {

			listCursor.fwd();

			Point cord = new Point(listCursor);

			if (listCursor.getDoublePosition(direction) < point[direction])

				Leftsublist.add(cord, listCursor.get());

			else

				Rightsublist.add(cord, listCursor.get());

		}
	
		

	final	Node<T> node = new Node<T>(point,direction,treeindex, Leftsublist, Rightsublist,null,null);
	
		if(node.Leftsublist.realMax(otherdirection)-node.Leftsublist.realMin(otherdirection)+1>2)
		node.leftchild = makeNode(node.Leftsublist, otherdirection, treeindex+1);
		if(node.Rightsublist.realMax(otherdirection)-node.Leftsublist.realMin(otherdirection)+1>2)
		node.rightchild = makeNode(node.Rightsublist,otherdirection, treeindex+1);
		

		return node;
	
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

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());

		int n = bitimg.numDimensions();

		FloatType val = new FloatType(200);

		// ImageJFunctions.show(img).setTitle("KD-Tree input");

		createBitimage(img, bitimg, val);

		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		IterableInterval<BitType> view = Views.interval(bitimg, new long[] { 0, 0 }, new long[] { 100, 100 });

		list = getList(view);

		PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(n);

		PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(n);

		listonlyones = getvalueList(list, 1);
		listonlyzeros = getvalueList(list, 0);

		

		
Node<BitType> testnode = makeNode(list, 0, 0);

Node<BitType> lefttestnode = testnode.leftchild;

System.out.print(lefttestnode.leftchild.treeindex);
	}

}
