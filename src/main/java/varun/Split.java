package varun;

import java.util.ArrayList;

import net.imglib2.Cursor;

import net.imglib2.IterableRealInterval;

import net.imglib2.RealCursor;

import net.imglib2.type.numeric.RealType;

public class Split {

	public static <T extends RealType<T>> void makeTree(final IterableRealInterval<T> interval) {
		int n = interval.numDimensions();
		long size = interval.size();
		double min[] = new double[n];
		double max[] = new double[n];

		final ArrayList<Cursor<T>> list = new ArrayList<Cursor<T>>((int) size);

		final RealCursor<T> cursor = interval.localizingCursor();

		while (cursor.hasNext()) {
			cursor.next();

		}

	}

}
