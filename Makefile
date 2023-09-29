feyfey:
	g++ -g main.cpp -o feyfey

feyfey-opti:
	g++ -g main.cpp -O3 -o feyfey

clean:
	rm -f feyfey