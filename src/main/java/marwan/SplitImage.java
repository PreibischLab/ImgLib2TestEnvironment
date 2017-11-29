package marwan;

import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class SplitImage {

	public <T extends NumericType<T> & NativeType<T>> SplitImage() {
		File file = new File("src/main/resources/DrosophilaWing.tif");
		final ImagePlus imp = new Opener().openImage(file.getAbsolutePath());

		// We will split the image in 4 columns 2 rows
		final int columns = 4;
		final int rows = 2;

		ArrayList<Rectangle> elements = splitImage(imp, columns, rows);
	}

	<T extends NumericType<T> & NativeType<T>> void openImage() {

	}

	ArrayList<Rectangle> splitImage(ImagePlus imp, int columns, int rows) {

		int totalWidth = imp.getWidth();
		int totalHeight = imp.getHeight();
		int width = Math.round(totalWidth / columns) - 1;
		int height = Math.round(totalHeight / rows) - 1;
		double overlap = 0.02;

		ArrayList<Rectangle> elements = new ArrayList<Rectangle>();

		// wrap it into an ImgLib image (no copying)
		final Img<FloatType> img = ImagePlusAdapter.wrap(imp);

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				Rectangle rect = new Rectangle(j * width, i * height, width, height);
				elements.add(rect);
				RandomAccessibleInterval<FloatType> view = Views.interval(img, new long[] { j * width, i * height },
						new long[] { (j + 1) * width, (i + 1) * height });
				ImageJFunctions.show(view);
			}
		}
		return elements;
	}

	public static void main(String[] args) {

		new ImageJ();
		new SplitImage();
	}
}
