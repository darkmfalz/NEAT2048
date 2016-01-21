import java.util.Random;

public class Game {
	
	private int[][] board;
	public int score = 0;
	public boolean canWin;
	
	public Game(){
		
		board = new int[4][4];
		canWin = true;
		
		addNum();
		addNum();
		
	}
	
	public int[][] getBoard(){
		
		return board;
		
	}
	
	public void resetBoard(){
		
		board = new int[4][4];
		score = 0;
		canWin = true;
		
		addNum();
		addNum();
		
	}
	
	private void addNum(){
		
		int a = 1;
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 0; j < board[i].length; j++){
				
				a = a*board[i][j];
				
				if(a >= 1)
					a = 1;
				
			}
			
		}
		
		if(a == 0){
			
			Random random = new Random();
		
			int x = random.nextInt(4);
			int y = random.nextInt(4);
			
			while(board[x][y] != 0){
				
				x = random.nextInt(4);
				y = random.nextInt(4);
				
			}
			
			boolean num = random.nextBoolean();
			
			if(num)
				board[x][y] = 2;
			else
				board[x][y] = 4;
		
		}
		else
			canWin = false;
			
	}
	
	public void printBoard(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 0; j < board[i].length; j++){
				
				System.out.print(board[i][j] + " ");
				
			}
			
			System.out.println();
			
		}
		
	}
	
	public int findMax(){
		
		int max = 0;
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 0; j < board[i].length; j++){
				
				if(board[i][j] > max)
					max = board[i][j];
				
			}
			
		}
		
		return max;
		
	}
	
	public int sumAll(){
		
		int sum = 0;
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 0; j < board[i].length; j++){
				
				sum += board[i][j];
				
			}
			
		}
		
		return sum;
		
	}
	
	public void move(int dir){
		
		int[][] boardCopy = board.clone();
		
		switch(dir){
			case 0:
				for(int i = 0; i < 4; i++)
					left();
				combineLeft();
				for(int i = 0; i < 4; i++)
					left();
				break;
			case 1:
				for(int i = 0; i < 4; i++)
					right();
				combineRight();
				for(int i = 0; i < 4; i++)
					right();
				break;
			case 2:
				for(int i = 0; i < 4; i++)
					up();
				combineUp();
				for(int i = 0; i < 4; i++)
					up();
				break;
			case 3:
				for(int i = 0; i < 4; i++)
					down();
				combineDown();
				for(int i = 0; i < 4; i++)
					down();
				break;
		}
		
		if(!board.equals(boardCopy))
			addNum();
		
	}
	
	public void randomMove(){
		
		Random random = new Random();
		move(random.nextInt(4));
		
	}
	
	public void playGameRando(int moves){
		
		for(int i = 0; i < moves; i++){
			
			if(canWin){
				
				randomMove();
				
			}
			else
				break;
				
		}
		
		resetBoard();
		
	}
	
	private void left(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = board[i].length - 2; j >= 0; j--){
				
				if(board[i][j + 1] > 0){
					
					if(board[i][j] == 0){
						
						board[i][j] = board[i][j + 1];
						board[i][j + 1] = 0;
						
					}
					
				}
				
			}
			
		}
		
	}

	private void combineLeft(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = board[i].length - 2; j >= 0; j--){
				
				if(board[i][j + 1] > 0){
					
					if(board[i][j] == board[i][j + 1]){
						
						board[i][j] = board[i][j] + board[i][j + 1]; 
						board[i][j + 1] = 0;
						
						score += board[i][j];
						
					}
					
				}
				
			}
			
		}
		
	}
	
	private void right(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 1; j < board[i].length; j++){
				
				if(board[i][j - 1] > 0){
					
					if(board[i][j] == 0){
						
						board[i][j] = board[i][j - 1];
						board[i][j - 1] = 0;
						
					}
					
				}
				
			}
			
		}
		
	}
	
	private void combineRight(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 1; j < board[i].length; j++){
				
				if(board[i][j - 1] > 0){
					
					if(board[i][j] == board[i][j - 1]){
						
						board[i][j] = board[i][j] + board[i][j - 1]; 
						board[i][j - 1] = 0;
						
						score += board[i][j];
						
					}
					
				}
				
			}
			
		}
		
	}
	
	private void up(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = board[i].length - 2; j >= 0; j--){
				
				if(board[j + 1][i] > 0){
					
					if(board[j][i] == 0){
						
						board[j][i] = board[j + 1][i];
						board[j + 1][i] = 0;
						
					}
					
				}
				
			}
			
		}
		
	}
	
	private void combineUp(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = board[i].length - 2; j >= 0; j--){
				
				if(board[j + 1][i] > 0){
					
					if(board[j][i] == board[j + 1][i]){
						
						board[j][i] = board[j][i] + board[j + 1][i]; 
						board[j + 1][i] = 0;
						
						score += board[j][i];
						
					}
					
				}
				
			}
			
		}
		
	}
	
	private void down(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 1; j < board[i].length; j++){
				
				if(board[j - 1][i] > 0){
					
					if(board[j][i] == 0){
						
						board[j][i] = board[j - 1][i];
						board[j - 1][i] = 0;
						
					}
					
				}
				
			}
			
		}
		
	}
	
	private void combineDown(){
		
		for(int i = 0; i < board.length; i++){
			
			for(int j = 1; j < board[i].length; j++){
				
				if(board[j - 1][i] > 0){
					
					if(board[j][i] == board[j - 1][i]){
						
						board[j][i] = board[j][i] + board[j - 1][i]; 
						board[j - 1][i] = 0;
						
						score += board[j][i];
						
					}
					
				}
				
			}
			
		}
		
	}

}
