package marwan;

import java.util.ArrayList;
import java.util.concurrent.Callable;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import marwan.Helper.Task;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class MultithreadingConvolution {

	public MultithreadingConvolution() throws ImgIOException, InterruptedException, IncompatibleTypeException {
		Helper.count = 0;
		Helper.sigma = 8;
		Helper.log = false;

		String string = "src/main/resources/mri-stack.tif";
//		String string = "src/main/resources/DrosophilaWing.tif";
		Img<FloatType> image = new ImgOpener().openImg(string, new FloatType());
	

		Img<FloatType> resultImage = ArrayImgs.floats(Helper.getDimensions(image));

		ImageJFunctions.show(image).setTitle("Origin");

		ImageJFunctions.show(resultImage).setTitle("final Initial");
		
		ArrayList<Portion> portions = new ArrayList<Portion>();
		Helper.splitImage(image, portions,-1,0);

		final ArrayList<myTask> taskList = Helper.createThreadTasks(portions, resultImage, Task.Gaus);
		for (myTask task : taskList)
			task.run();

		ImageJFunctions.show(resultImage).setTitle("final");
	}

	public static void main(String[] args) throws ImgIOException, InterruptedException, IncompatibleTypeException {

		new ImageJ();
		new MultithreadingConvolution();
	}
}
