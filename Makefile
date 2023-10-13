feyfey:
	g++ -g main.cpp -o feyfey

feyfey-opti:
	g++ -g main.cpp -O3 -o feyfey

divideandconquer:
	g++ -g divideandconquer.cpp -o divideandconquer

divideandconquer-opti:
	g++ -g divideandconquer.cpp -O3 -o divideandconquer

clean:
	rm -f feyfey divideandconquer