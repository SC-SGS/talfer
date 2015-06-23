#ifndef GAME_VERSION2_HPP_
#define GAME_VERSION2_HPP_

#include <time.h>
#include <sstream>
#include <string>
#include "SDL.h"
#include "SDL_mixer.h"
#include "CFont.hpp"
#include "CStage_VideoOutput.hpp"
#include "CParticle.hpp"

class Game {
private:
	int points;
	time_t start_time, end_time;
	int elapse_seconds;
	/* Mix_Music actually holds the music information.  */

	Mix_Music *music;
	Mix_Chunk *bumped;
	Mix_Chunk *destroyed;
	Mix_Chunk *points_sound;
	Mix_Chunk *spawn_sound;

	CStage_VideoOutput* vo;
	std::string points_s;

	static const int small_penguin = 6;
	static const int big_penguin = 14;

public:
	bool global_pause;
	bool has_ended;

public:
	static Game* getInstance() {
		static Game instance;

		return &instance;
	}

private:
	Game() :
			points(0), bumped(0), destroyed(0), points_sound(0), spawn_sound(0), points_s(
					"Points: 0"), global_pause(false), has_ended(false) {
		// Get game start time
		time(&start_time);
		init_audio();
		load_sounds();
		play_music();
	}

private:
	// Do not implement:
	Game(Game const&);
	void operator=(Game const&);

private:
	void load_sounds() {
		bumped = Mix_LoadWAV("data/music/collision.wav");
		if (!bumped) {
			printf("Mix_LoadWAV(\"collision.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		} else {
			Mix_VolumeChunk(bumped, 32);
		}

		destroyed = Mix_LoadWAV("data/music/destroyed.wav");
		if (!destroyed) {
			printf("Mix_LoadWAV(\"destroyed.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		}

		points_sound = Mix_LoadWAV("data/music/points.wav");
		if (!points_sound) {
			printf("Mix_LoadWAV(\"points.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		}

		spawn_sound = Mix_LoadWAV("data/music/spawn.wav");
		if (!spawn_sound) {
			printf("Mix_LoadWAV(\"spawn.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		} else {
			Mix_VolumeChunk(spawn_sound, 32);
		}
	}

public:
	void particle_destroyed(float mass, int bumps, float alive_time,
			bool target, int x, int y) {
		int new_points = mass / (bumps + 1) * alive_time;

		if (!(target && mass == small_penguin * small_penguin)
				|| !(!target && mass == big_penguin * big_penguin)) {
			// Giving points
			points += new_points;

			// TODO Give feedback on amount of points
			if (Mix_PlayChannel(3, destroyed, 0) == -1) {
				printf("Mix_PlayChannel: %s\n", Mix_GetError());
				// well, there's no music, but most games don't break without music...
			}

		} else if (!(target && mass != small_penguin * small_penguin)
				|| !(!target && mass == small_penguin * small_penguin)) {

			// Giving points
			new_points = -new_points;
			points += new_points;

			// TODO Give feedback on amount of points
			if (Mix_PlayChannel(3, destroyed, 0) == -1) {
				printf("Mix_PlayChannel: %s\n", Mix_GetError());
				// well, there's no music, but most games don't break without music...
			}
		}

		// Update point string
		std::stringstream ss;
		ss.flush();
		ss << "Points: " << points;
		points_s = ss.str();
		vo->change_Points(&points_s);

		// Show +Points
		std::stringstream ss2;
		ss2.flush();
		if (new_points < 0) {
			ss2 << (int) new_points;
		} else if (new_points >= 0) {
			ss2 << "+" << (int) new_points;
		}
		vo->push_Text(ss2.str(), 15, x, y, 20);

		/*			if (collision != nullptr && particles.at(i)->getNextParticleRadius()==6){
		 float lifetime = (clock() - particles.at(i)->spawn_time) / CLOCKS_PER_SEC;
		 game->particle_destroyed(1 / (float) particles.at(i)->inverseMass, lifetime, false);
		 particles.erase(particles.begin()+i);
		 i--;
		 }
		 else if (collision != nullptr && particles.at(i)->getNextParticleRadius()==3){
		 float lifetime = (clock() - particles.at(i)->spawn_time) / CLOCKS_PER_SEC;
		 game->particle_destroyed(1 / (float) particles.at(i)->inverseMass, lifetime, true);
		 particles.erase(particles.begin()+i);
		 i--;
		 }
		 else{
		 collision = checkFlagArea(particles.at(i)->position, 4, particles.at(i)->radius);
		 }

		 if (collision != nullptr && particles.at(i)->getNextParticleRadius()==6){
		 float lifetime = (clock() - particles.at(i)->spawn_time) / CLOCKS_PER_SEC;
		 game->particle_destroyed(1 / (float) particles.at(i)->inverseMass, lifetime, true);
		 particles.erase(particles.begin()+i);
		 i--;
		 }
		 else if (collision != nullptr && particles.at(i)->getNextParticleRadius()==3){
		 float lifetime = (clock() - particles.at(i)->spawn_time) / CLOCKS_PER_SEC;
		 game->particle_destroyed(1 / (float) particles.at(i)->inverseMass, lifetime, false);
		 particles.erase(particles.begin()+i);
		 i--;
		 }

		 // Calc new points to be awarded
		 // target: "big particle" in blue outflow or "small particle" in green outflow
		 // -> giving points
		 int new_points;
		 if (target){
		 new_points = mass/(1)*alive_time;
		 }
		 // target: "big particle" in green outflow or "small particle" in blue outflow
		 // -> giving negativ points
		 else if (!target){
		 new_points = (-1)*(mass/(1)*alive_time)/3;
		 }
		 // target: particle doesn't reach an outflow
		 // -> giving no points
		 else {
		 new_points = 0;
		 }

		 // Giving points
		 points += new_points;

		 // TODO Give feedback on amount of points
		 if(Mix_PlayChannel(3, destroyed, 0) == -1) {
		 printf("Mix_PlayChannel: %s\n", Mix_GetError());
		 // well, there's no music, but most games don't break without music...
		 }

		 // Update point string
		 std::cout << "Points: " << points << "\t+" << mass/(1)*alive_time << std::endl;

		 */}

	void particle_bumped() {
		if (Mix_PlayChannel(2, bumped, 0) == -1) {
			printf("Mix_PlayChannel: %s\n", Mix_GetError());
			// well, there's no music, but most games don't break without music...
		}
	}

	// particle mass is then calculated as radius^2
	void setNextParticle(CParticle * p) {
		if (rand() % 2 == 1) {
			p->radius = small_penguin;
			p->color = 1;
		} else {
			p->radius = big_penguin;
			p->color = 0;
		}
	}

	void speed(float min, float max) {
		// TODO
		// give feedback
		std::cout << "Particle Speeds min:" << min << " max:" << max
				<< std::endl;
	}

	void gesture(std::string /*name*/, float /*direction_x*/,
			float /*direction_y*/) {
		// TODO
		// Displays gesture recognition feedback
	}

	void particle_created(float /*mass*/) {
		if (Mix_PlayChannel(1, spawn_sound, 0) == -1) {
			printf("Mix_PlayChannel: %s\n", Mix_GetError());
			// well, there's no music, but most games don't break without music...
		}
	}

	static void end_game() {
		Game* g = getInstance();

		// Get game end time
		time(&g->end_time);
		// Calculate elapsed seconds
		g->elapse_seconds = difftime(g->start_time, g->end_time);

		g->vo->push_Text("Game Over!", -1, 75, 250, 100);

		// Halt game
		g->global_pause = true;

		g->has_ended = true;
	}

private:
	void init_audio() {
		/* We're going to be requesting certain things from our audio
		 device, so we set them up beforehand */
		int audio_rate = 22050;
		Uint16 audio_format = AUDIO_S16; /* 16-bit stereo */
		int audio_channels = 2;
		int audio_buffers = 4096;

		/* This is where we open up our audio device.  Mix_OpenAudio takes
		 as its parameters the audio format we'd /like/ to have. */
		if (Mix_OpenAudio(audio_rate, audio_format, audio_channels,
				audio_buffers)) {
			printf("Unable to open audio!\n");
			exit(1);
		}

		// allocate 16 mixing channels
		Mix_AllocateChannels(16);
	}

	void play_music() {
		/* Actually loads up the music */
		music = Mix_LoadMUS("data/music/bg-music.wav");
		if (!music) {
			printf("Mix_LoadMUS(\"bg-music.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		}

		/* This begins playing the music - the first argument is a
		 pointer to Mix_Music structure, and the second is how many
		 times you want it to loop (use -1 for infinite, and 0 to
		 have it just play once) */
		if (Mix_PlayMusic(music, 0) == -1) {
			printf("Mix_PlayMusic: %s\n", Mix_GetError());
			// well, there's no music, but most games don't break without music...
		}
		Mix_HookMusicFinished(end_game);
	}
public:
	void setVideoOutput(CStage_VideoOutput* vo) {
		this->vo = vo;

		vo->change_Points(&points_s);
		//std::cout<< points_s << std::endl;
	}
public:
	void mute_music(bool mute) {
		if (mute) {
			// Turn down volume
			Mix_VolumeMusic(0);
		} else {
			// Turn up volume
			Mix_VolumeMusic(128);
		}
	}
};

#endif
