public class Main {

	public static void main(String[] args){
		
		Game game = new Game();
		game.printBoard();
		System.out.println();
		
		for(int i = 0; i < 100; i++)
			game.playGameRando(100);
		
		Breeder breeder = new Breeder(100);
		
	}
	
}
