# Documentation API DuoBingo

## Vue d'ensemble

L'API DuoBingo fournit un ensemble d'endpoints RESTful pour gérer l'apprentissage de langues, les utilisateurs, les exercices et le suivi de progression.

## URL de base

```
https://api.duobingo.com/v1
```

## Authentification

L'API utilise l'authentification par token JWT (JSON Web Token). Pour accéder aux endpoints protégés, incluez le token dans l'en-tête de la requête :

```
Authorization: Bearer <votre_token_jwt>
```

### Obtenir un token

**Endpoint:** `POST /auth/login`

**Corps de la requête:**
```json
{
  "email": "utilisateur@example.com",
  "password": "motdepasse"
}
```

**Réponse:**
```json
{
  "token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...",
  "user": {
    "id": "123",
    "email": "utilisateur@example.com",
    "username": "utilisateur"
  }
}
```

## Endpoints principaux

### Utilisateurs

#### Créer un utilisateur

- **Endpoint:** `POST /users`
- **Authentification:** Non requise
- **Corps de la requête:**

```json
{
  "username": "nouveau_utilisateur",
  "email": "user@example.com",
  "password": "motdepasse_securise"
}
```

- **Réponse (201):**

```json
{
  "id": "456",
  "username": "nouveau_utilisateur",
  "email": "user@example.com",
  "createdAt": "2024-01-14T10:00:00Z"
}
```

#### Obtenir un utilisateur

- **Endpoint:** `GET /users/{userId}`
- **Authentification:** Requise
- **Réponse (200):**

```json
{
  "id": "456",
  "username": "nouveau_utilisateur",
  "email": "user@example.com",
  "level": 5,
  "points": 1250,
  "createdAt": "2024-01-14T10:00:00Z"
}
```

### Leçons

#### Lister les leçons

- **Endpoint:** `GET /lessons`
- **Authentification:** Requise
- **Paramètres de requête:**
  - `language` (optionnel): Code de langue (ex: "en", "fr", "es")
  - `level` (optionnel): Niveau de difficulté (1-10)

- **Réponse (200):**

```json
{
  "lessons": [
    {
      "id": "lesson_1",
      "title": "Les bases de la conversation",
      "language": "fr",
      "level": 1,
      "description": "Apprenez les phrases de base pour commencer une conversation"
    }
  ]
}
```

#### Obtenir une leçon

- **Endpoint:** `GET /lessons/{lessonId}`
- **Authentification:** Requise
- **Réponse (200):**

```json
{
  "id": "lesson_1",
  "title": "Les bases de la conversation",
  "language": "fr",
  "level": 1,
  "description": "Apprenez les phrases de base pour commencer une conversation",
  "exercises": [
    {
      "id": "ex_1",
      "type": "translation",
      "question": "Bonjour",
      "answer": "Hello"
    }
  ]
}
```

### Exercices

#### Soumettre une réponse

- **Endpoint:** `POST /exercises/{exerciseId}/submit`
- **Authentification:** Requise
- **Corps de la requête:**

```json
{
  "answer": "Hello"
}
```

- **Réponse (200):**

```json
{
  "correct": true,
  "points": 10,
  "feedback": "Excellente réponse !",
  "correctAnswer": "Hello"
}
```

### Progression

#### Obtenir la progression de l'utilisateur

- **Endpoint:** `GET /users/{userId}/progress`
- **Authentification:** Requise
- **Réponse (200):**

```json
{
  "userId": "456",
  "totalPoints": 1250,
  "level": 5,
  "completedLessons": 15,
  "streak": 7,
  "languages": [
    {
      "code": "fr",
      "name": "Français",
      "level": 5,
      "points": 1250
    }
  ]
}
```

## Codes d'erreur

L'API utilise les codes de statut HTTP standard :

- **200 OK** - Requête réussie
- **201 Created** - Ressource créée avec succès
- **400 Bad Request** - Requête invalide
- **401 Unauthorized** - Authentification requise ou invalide
- **403 Forbidden** - Accès interdit
- **404 Not Found** - Ressource non trouvée
- **500 Internal Server Error** - Erreur serveur

### Format d'erreur

```json
{
  "error": {
    "code": "INVALID_REQUEST",
    "message": "Description de l'erreur",
    "details": []
  }
}
```

## Limites de taux

- **Requêtes authentifiées:** 1000 requêtes par heure
- **Requêtes non authentifiées:** 100 requêtes par heure

## Pagination

Les endpoints qui retournent des listes supportent la pagination :

**Paramètres:**
- `page` (défaut: 1) - Numéro de page
- `limit` (défaut: 20, max: 100) - Nombre d'éléments par page

**Exemple:**
```
GET /lessons?page=2&limit=10
```

**Réponse:**
```json
{
  "data": [...],
  "pagination": {
    "page": 2,
    "limit": 10,
    "total": 45,
    "pages": 5
  }
}
```

## Support

Pour plus d'informations, consultez la [spécification OpenAPI](./openapi.yaml) complète.
