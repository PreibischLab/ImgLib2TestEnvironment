package marwan;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class MultithreadingConvolution {
	

public MultithreadingConvolution() throws ImgIOException, InterruptedException {
	Helper.count = 0;
	Helper.log = false;
		String string = "src/main/resources/DrosophilaWing.tif";
		@SuppressWarnings("deprecation")
		Img<FloatType> image = new ImgOpener().openImg(string, new FloatType());
	
		ImageJFunctions.show(image);
//		Img< FloatType > resultImage = ArrayImgs.floats( 1024, 1024 );
//		TODO change it to variable
//		int numTasks = 20;
		
		
		final ArrayList<myTask> taskList = createThreadTasks(image);
		for(myTask task:taskList) task.run();
		
		ImageJFunctions.show(image);
		
	}


	
	private ArrayList<myTask> createThreadTasks(RandomAccessibleInterval<FloatType> image) {
		final int columns = 4;
		final int rows = 2;

		ArrayList<RandomAccessibleInterval<FloatType>> views = Helper.splitImage(image, columns, rows);
		ArrayList<myTask> taskList = new ArrayList<myTask>();
		RandomAccessible<FloatType> infiniteImg = Views.extendMirrorSingle(image);
		 for(RandomAccessibleInterval<FloatType> view: views) {
			 myTask task = new myTask(infiniteImg, view);
			 taskList.add(task);
		 }
		
	return taskList;
}



	ArrayList< Callable< float[] > > createThreadTasks(Img<FloatType> image, Img<FloatType> resultImage, int numTasks, int stepSize) {
		ArrayList< Callable< float[] > > taskList = new ArrayList< Callable< float[] > >();
		for (int i = 0; i < numTasks; i++) {

			int startPosition = i*stepSize;
			Callable<float[]> callable = new Callable<float[]>() {

				@Override
				public float[] call() throws Exception {
					float[] result = traitment(image,resultImage,numTasks,stepSize);
					return result;
				}
				
			};
			taskList.add(callable);
			
		}
		return taskList;
		
	}
	

	protected float[] traitment(Img<FloatType> image, Img<FloatType> resultImage, int numTask, int stepSize) {
	 float min = Helper.computeMiLocation(image, numTask*stepSize, (numTask+1)*stepSize);
	return new float[] {min};
		
	}

	public static void main(String[] args) throws ImgIOException, InterruptedException {

		new ImageJ();
		new MultithreadingConvolution();
	}
}
