
package varun;
import java.io.File;
import java.io.FileNotFoundException;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;
import varun.MyKDtree.EucledianDistance;

public class SortList{
	public static <T extends RealType<T>> void split(PointSampleList<T> list, int direction) {

		int n = list.numDimensions();

		if (list.size() <= 1)
			return;

		else {
/****
			 * To ward against running over the dimensionality, creating some
			 * local restrictions on the global variable direction
			 ****/
			if (direction == list.numDimensions())
				direction = 0;

			// the first element belonging to the right list childB
			final int splitIndex = (int) list.size() / 2;

			final PointSampleList<T> childA = new PointSampleList<T>(n);
			final PointSampleList<T> childB = new PointSampleList<T>(n);

			final Cursor<T> listCursor = list.localizingCursor();

			int index = 0;
			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(listCursor);

				if (index < splitIndex) {

					childA.add(cord, listCursor.get().copy());

					// System.out.println("childA: "+listCursor.get());

				} else

				{

					childB.add(cord, listCursor.get().copy());
					// System.out.println("childB: "+listCursor.get());
				}
				index++;
			}

			split(childA, direction);

			split(childB, direction);

			mergeList(list, childA, childB, direction);
		}

	}

	/// ***** Returns a sorted list *********////
	public static <T extends RealType<T>> void mergeList(PointSampleList<T> list, PointSampleList<T> listA,
			PointSampleList<T> listB, int direction) {

		final Cursor<T> cursorA = listA.localizingCursor();
		final Cursor<T> cursorB = listB.localizingCursor();
		final Cursor<T> cursor = list.localizingCursor();

		cursorA.fwd();
		cursorB.fwd();

		boolean cannotMoveOn = false;

		do {
			// here is where you decide what you sort after
			if (cursorA.getDoublePosition(direction) < (cursorB.getDoublePosition(direction))) {

				cursor.fwd();
				cursor.get().set(cursorA.get());
				if (cursorA.hasNext())
					cursorA.fwd();
				else {
					cannotMoveOn = true;

					// move cursorB until the end
					boolean stopped = false;
					do {
						cursor.fwd();
						cursor.get().set(cursorB.get());
						if (cursorB.hasNext())
							cursorB.fwd();
						else
							stopped = true;
					} while (stopped == false);
				}

			}

			else

			{

				cursor.fwd();
				cursor.get().set(cursorB.get());
				if (cursorB.hasNext())
					cursorB.fwd();
				else {
					cannotMoveOn = true;

					// move cursorA until the end
					boolean stopped = false;
					do {
						cursor.fwd();
						cursor.get().set(cursorA.get());
						if (cursorA.hasNext())
							cursorA.fwd();
						else
							stopped = true;
					} while (stopped == false);
				}

			}

		} while (cannotMoveOn == false);
	}

	
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


	public static void main(String[] args) throws FileNotFoundException {
		 long startTime = System.currentTimeMillis();
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(img, new BitType());
		final Img<FloatType> imgout = new ArrayImgFactory<FloatType>().create(img, new FloatType());
		PointSampleList<BitType> list = new PointSampleList<BitType>(bitimg.numDimensions());

		IterableInterval<BitType> view = Views.interval(bitimg, new long[] { 0, 0 }, new long[] { 127, 127 });

		list = getList(view);

		PointSampleList<BitType> listonlyones = new PointSampleList<BitType>(bitimg.numDimensions());

		PointSampleList<BitType> listonlyzeros = new PointSampleList<BitType>(bitimg.numDimensions());

		listonlyones = getvalueList(list, 1);
		listonlyzeros = getvalueList(list, 0);

		split(listonlyones, 0);
		split(listonlyones, 1);

	
		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println(totalTime);
		
		
	}
}
