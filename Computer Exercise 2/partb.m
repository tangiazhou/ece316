% Computer Exercise 2

% Part B
% Modify the above program to create a Morse code signal which
% transmits: Hello Elvino is Here.
% Substitute "Elvino" with your first name.
% Look up the morse code here
% https://en.wikipedia.org/wiki/Morse_code#/media/File:International_Morse_Code.svg
% Note: Each letter is represented by a sequence of "dots" and "dashes". The length of the 
% dash equals 3 times the lenght of the dot. Use a length of one dot as the space between dots
% and dashes (i.e. quiet). Use the length of 4 dots as the space between two letters.
% Choose a suitable length of time for the "dot".

clear all;
close all;

% Morse code program that plays "Hello Tangia is Here"

t0 = 0:5000; % Define a time vector

dot = cos(t0(1:1000));
dash = cos(t0(1:3000)); % Dash 3 times the length of the dot 
silence1 = zeros(1,1000); % Silence between dots and dashes
silence2 = zeros(1, 4000); % Silence between letters 

h = [dot silence1 dot silence1 dot silence1 dot silence2];
e = [dot silence2];
l = [dot silence1 dash silence1 dot silence1 dot silence2];
o = [dash silence1 dash silence1 dash silence2];
t = [dash silence2];
a = [dot silence1 dash silence2];
n = [dash silence1 dot silence2];
g = [dash silence1 dash silence1 dot silence2];
i = [dot silence1 dot silence2];
s = [dot silence1 dot silence1 dot silence2];
r = [dot silence1 dash silence1 dot silence2];

% H
player =audioplayer(h,11000);
play(player);

% E
player=audioplayer(e,11000);
play(player);

% L
player=audioplayer(l,11000);
play(player);

% L
player=audioplayer(l,11000);
play(player);

% O
player=audioplayer(o,11000);
play(player);

% T
player=audioplayer(t,11000);
play(player);

% A
player=audioplayer(a,11000);
play(player);

% N
player=audioplayer(n,11000);
play(player);

% G
player=audioplayer(g,11000);
play(player);

% I
player=audioplayer(i,11000);
play(player);

% A
player=audioplayer(a,11000);
play(player);

% I
player=audioplayer(i,11000);
play(player);

% S
player=audioplayer(s,11000);
play(player);

% H
player=audioplayer(h,11000);
play(player);

% E
player=audioplayer(e,11000);
play(player);

% R
player=audioplayer(r,11000);
play(player);

% E
player=audioplayer(e,11000);
play(player);
