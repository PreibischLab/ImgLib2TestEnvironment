package marwan;

import java.util.ArrayList;
import java.util.concurrent.Callable;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import marwan.Helper.Task;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
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
		final int columns = 10;
		final int rows = 10;
		
		
		String string = "src/main/resources/DrosophilaWing.tif";
		Img<FloatType> image = new ImgOpener().openImg(string, new FloatType());

		Img<FloatType> resultImage = ArrayImgs.floats(Helper.getDimensions(image));
		RandomAccessible<FloatType> infiniteImg = Views.extendMirrorSingle(image);

	
		FinalInterval interval = Helper.getFinalInterval(image);
		ImageJFunctions.show(image).setTitle("Origin");
		ImageJFunctions.show(Views.interval(infiniteImg, interval)).setTitle("Extented Image");

		ImageJFunctions.show(resultImage).setTitle("final Initial");

	
		ArrayList<Portion> portions = Helper.splitImage(infiniteImg,interval, columns, rows);
		final ArrayList<myTask> taskList = createThreadTasks(portions, resultImage,Task.Gaus);
		for (myTask task : taskList)
			task.run();


ImageJFunctions.show(resultImage).setTitle("final");
	}


	ArrayList<myTask> createThreadTasks(ArrayList<Portion> portions,Img<FloatType> resultImage, Task type) {
		
		ArrayList<myTask> taskList = new ArrayList<myTask>();
		switch (type) {
		case Gaus:
			for (Portion portion : portions) {
				myTask task = new myTask(portion,resultImage,Task.Gaus);
				taskList.add(task);
			}
			break;

		default:
			break;
		}
		return taskList;
	}

	ArrayList<Callable<float[]>> createThreadTasks(Img<FloatType> image, Img<FloatType> resultImage, int numTasks,
			int stepSize) {
		ArrayList<Callable<float[]>> taskList = new ArrayList<Callable<float[]>>();
		for (int i = 0; i < numTasks; i++) {

			int startPosition = i * stepSize;
			Callable<float[]> callable = new Callable<float[]>() {

				@Override
				public float[] call() throws Exception {
					float[] result = traitment(image, resultImage, numTasks, stepSize);
					return result;
				}

			};
			taskList.add(callable);

		}
		return taskList;

	}

	protected float[] traitment(Img<FloatType> image, Img<FloatType> resultImage, int numTask, int stepSize) {
		float min = Helper.computeMiLocation(image, numTask * stepSize, (numTask + 1) * stepSize);
		return new float[] { min };

	}

	public static void main(String[] args) throws ImgIOException, InterruptedException, IncompatibleTypeException {

		new ImageJ();
		new MultithreadingConvolution();
	}
}
