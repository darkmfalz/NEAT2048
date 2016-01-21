public class Main {

	public static void main(String[] args){
		
		BreederAlt breeder = new BreederAlt(1000);

		while(breeder.maxFitness < 2500)
			breeder.breedNextGen();
		
	}
	
}
