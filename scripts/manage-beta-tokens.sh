#!/bin/bash
# CABS-flex Beta Token Management Script
# Generate and manage tokens for beta testers

set -e

TOKENS_FILE="beta-tokens.txt"
TOKEN_PREFIX="cabsflex_beta_"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

show_usage() {
    echo "CABS-flex Beta Token Manager"
    echo "Usage: $0 [command] [options]"
    echo ""
    echo "Commands:"
    echo "  generate <email>     Generate new beta token"
    echo "  list                 List all beta tokens"
    echo "  revoke <token>       Revoke a beta token"
    echo "  instructions <token> Show installation instructions"
    echo ""
}

generate_token() {
    local email=$1
    if [ -z "$email" ]; then
        echo "Error: Email required"
        echo "Usage: $0 generate user@example.com"
        exit 1
    fi
    
    # Generate random token
    local token="${TOKEN_PREFIX}$(openssl rand -hex 16)"
    local date=$(date '+%Y-%m-%d %H:%M:%S')
    
    # Save to file
    echo "$date,$email,$token,active" >> $TOKENS_FILE
    
    echo -e "${GREEN}✅ Beta token generated for $email${NC}"
    echo -e "${BLUE}Token: $token${NC}"
    echo ""
    echo "Installation instructions:"
    echo "curl -sSL https://raw.githubusercontent.com/LCBio/cabsflex/main/install-beta.sh | bash -s $token"
}

list_tokens() {
    if [ ! -f "$TOKENS_FILE" ]; then
        echo "No tokens found."
        return
    fi
    
    echo -e "${BLUE}Beta Tokens:${NC}"
    echo "Date                 Email                    Token                                      Status"
    echo "-------------------- ------------------------ ------------------------------------------ --------"
    
    while IFS=',' read -r date email token status; do
        printf "%-20s %-24s %-42s %s\n" "$date" "$email" "$token" "$status"
    done < "$TOKENS_FILE"
}

revoke_token() {
    local token=$1
    if [ -z "$token" ]; then
        echo "Error: Token required"
        echo "Usage: $0 revoke <token>"
        exit 1
    fi
    
    if [ ! -f "$TOKENS_FILE" ]; then
        echo "No tokens found."
        return
    fi
    
    # Update token status to revoked
    sed -i "s/,$token,active/,$token,revoked/g" "$TOKENS_FILE"
    echo -e "${YELLOW}Token revoked: $token${NC}"
}

show_instructions() {
    local token=$1
    if [ -z "$token" ]; then
        echo "Error: Token required"
        echo "Usage: $0 instructions <token>"
        exit 1
    fi
    
    echo -e "${BLUE}CABS-flex Beta Installation Instructions${NC}"
    echo "========================================"
    echo ""
    echo "Send this to your beta tester:"
    echo ""
    echo "---"
    echo "Thank you for participating in CABS-flex beta testing!"
    echo ""
    echo "To install CABS-flex Beta, run this single command:"
    echo ""
    echo "curl -sSL https://raw.githubusercontent.com/LCBio/cabsflex/main/install-beta.sh | bash -s $token"
    echo ""
    echo "Prerequisites:"
    echo "- Anaconda or Miniconda installed"
    echo "- Git installed"
    echo "- Linux or macOS (Windows via WSL)"
    echo ""
    echo "After installation, activate and test:"
    echo "conda activate cabs"
    echo "CABSflex --help"
    echo ""
    echo "Report issues: https://github.com/LCBio/cabsflex/issues"
    echo "Contact: k.wroblewski7@uw.edu.pl"
    echo "---"
}

# Main command processing
case $1 in
    generate)
        generate_token $2
        ;;
    list)
        list_tokens
        ;;
    revoke)
        revoke_token $2
        ;;
    instructions)
        show_instructions $2
        ;;
    *)
        show_usage
        exit 1
        ;;
esac
