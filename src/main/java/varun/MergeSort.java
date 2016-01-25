package varun;



import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;

public class MergeSort {
	public static void sort() {

	}

	public static int[] sort(int[] list) {

		
		int tmp;

		for (int i = 0; i < list.length; ++i){
		for (int j = i+1; j < list.length; ++j) {
			if (list[i] >= list[j]){
				tmp = list[i];
				list[i]=list[j];
			 list[j]=tmp;
		}

			
		}
		}

		return list;
	}

	public static int[] merge(int[] listA, int[] listB) {

		int [] list= new int[listA.length+listB.length];
		int j = 0;
		for (int i = 0; i < listA.length; ++i) {

			list[j] = listA[i];
			++j;
		}

		for (int i = 0; i < listB.length; ++i) {

			list[j] = listB[i];
++j;
		}

		return list;

	}

	public static Pair<int[], int[]> split(int[] list) {
		if (list.length <= 1)
			return null;

		int[] out1 = new int[list.length / 2];
		int[] out2 = new int[list.length / 2 + list.length % 2];

		int j = 0;
		for (int i = 0; i < out1.length; ++i) {
			out1[i] = list[j];
			++j;
		}

		for (int i = 0; i < out2.length; ++i) {
			out2[i] = list[j];
			++j;
		}

		Pair<int[], int[]> pair = new ValuePair<int[], int[]>(out1, out2);

		return pair;

	}

	public static int[] mergeSort(int[] list) {

		int[] finallist = new int[list.length];
		
		Pair<int[], int[]> pair;

		if (list.length > 1) {

			int[] listA = new int[list.length/2];
			int[] listB = new int[list.length/2+list.length%2];
			int[] sortedlist = new int[list.length];

			pair = split(list);
			listA = sort(pair.getA());
			listB = sort(pair.getB());

			

			sortedlist = merge(listA,listB);
			
			finallist=sort(sortedlist);

		}

		return finallist;

	}

	public static void main(String[] args) {

		int[] list = new int[8];
		int[] sortedlist = new int[8];

		list[0] = 5;
		list[1] = 2;
		list[2] = 4;
		list[3] = 7;
		list[4] = 1;
		list[5] = 3;
		list[6] = 2;
		list[7] = 6;

		sortedlist=mergeSort(list);
		
		for (int i=0; i<list.length; ++i){
			
			System.out.print(sortedlist[i]);
		}
		

	}

}
