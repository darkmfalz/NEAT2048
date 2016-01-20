public class Main {

	public static void main(String[] args){
		
		Breeder breeder = new Breeder(100);
		
		while(breeder.maxFitness < 2048)
			breeder.breedNextGen();
		
	}
	
}
