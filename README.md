# DuoBingo API

Bienvenue dans le projet DuoBingo API - une API RESTful pour gérer des sessions d'apprentissage de langues interactives.

## Description

DuoBingo est une plateforme d'apprentissage de langues qui permet aux utilisateurs de pratiquer et d'améliorer leurs compétences linguistiques à travers des exercices interactifs et des sessions de jeu.

## Fonctionnalités

- Gestion des utilisateurs et authentification
- Sessions d'apprentissage personnalisées
- Exercices interactifs de vocabulaire et grammaire
- Suivi de progression
- Système de points et récompenses

## Documentation

- [Documentation API complète](./API_DOC_DUOBINGO.md)
- [Spécification OpenAPI](./openapi.yaml)

## Démarrage rapide

### Prérequis

- Node.js >= 16.x
- npm >= 8.x

### Installation

```bash
npm install
```

### Configuration

Créez un fichier `.env` à la racine du projet avec les variables d'environnement nécessaires :

```
PORT=3000
DATABASE_URL=your_database_url
JWT_SECRET=your_secret_key
```

### Démarrage

```bash
npm start
```

L'API sera accessible sur `http://localhost:3000`

## Structure du projet

```
.
├── README.md
├── API_DOC_DUOBINGO.md
├── openapi.yaml
├── LICENSE
└── .gitignore
```

## Technologies utilisées

- Node.js
- Express.js
- OpenAPI 3.0

## Licence

Ce projet est sous licence MIT. Voir le fichier [LICENSE](./LICENSE) pour plus de détails.

## Contact

Pour toute question ou suggestion, veuillez ouvrir une issue sur ce dépôt.
